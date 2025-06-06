#!/usr/bin/env python
"""
Machine learning prediction of superconducting Tc using RBT features.

This script implements a graph neural network (GNN) that learns to predict
superconducting critical temperatures from RBT ledger features and structural
descriptors. The model can be trained on known superconductors and used to
rank new candidates.

Usage:
    # Train model
    python 02_ml_predict.py --train --supercon_data data/supercon.csv --rbt_features results/scored_candidates.csv

    # Predict Tc for new candidates  
    python 02_ml_predict.py --predict --model models/best_model.pt --candidates results/shortlist.csv --out results/ml_predictions.csv

Reference: Combines RBT theory with modern deep learning for materials discovery
"""

import argparse
import json
import pickle
from pathlib import Path
from typing import List, Dict, Any, Tuple, Optional
import numpy as np
import pandas as pd
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# ML libraries
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import Adam, AdamW
from torch.optim.lr_scheduler import ReduceLROnPlateau
from sklearn.model_selection import train_test_split, KFold
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.ensemble import RandomForestRegressor

# Graph ML
try:
    from torch_geometric.data import Data, DataLoader
    from torch_geometric.nn import GCNConv, GATConv, global_mean_pool, global_max_pool
    TORCH_GEOMETRIC_AVAILABLE = True
except ImportError:
    print("Warning: torch_geometric not available. Using fallback models.")
    TORCH_GEOMETRIC_AVAILABLE = False

# Materials science
from pymatgen.core import Composition, Structure
from pymatgen.analysis.local_env import CrystalNN
from dscribe.descriptors import SOAP, MBTR

from utils.ledger import rbt_superconductor_score, predict_tc_estimate
from utils.io import structure_from_composition, load_supercon_data


class RBTFeatureExtractor:
    """Extract RBT-specific features for ML models."""
    
    def __init__(self):
        self.scaler = StandardScaler()
        self.fitted = False
    
    def extract_composition_features(self, composition: Composition) -> np.ndarray:
        """Extract composition-based features."""
        features = []
        
        # Basic composition features
        features.extend([
            len(composition.elements),
            composition.num_atoms,
            composition.weight,
            composition.weight / composition.num_atoms,  # avg atomic mass
        ])
        
        # Electronic features
        metal_count = sum(1 for el in composition.elements if el.is_metal)
        features.extend([
            metal_count / len(composition.elements),  # metal fraction
            sum(el.X for el in composition.elements if el.X) / len(composition.elements),  # avg electronegativity
        ])
        
        # RBT-specific features would be added by the ledger scoring
        return np.array(features)
    
    def extract_structure_features(self, structure: Structure) -> np.ndarray:
        """Extract structure-based features using SOAP descriptors."""
        try:
            # SOAP descriptor for local environments
            soap = SOAP(
                species=list(set(site.specie.symbol for site in structure)),
                r_cut=6.0,
                n_max=8,
                l_max=6,
                periodic=True
            )
            
            # Get SOAP features for all sites
            soap_features = soap.create(structure)
            
            # Aggregate features (mean and std across sites)
            if len(soap_features.shape) > 1:
                features = np.concatenate([
                    np.mean(soap_features, axis=0),
                    np.std(soap_features, axis=0)
                ])
            else:
                features = soap_features
                
        except Exception as e:
            print(f"Warning: SOAP feature extraction failed: {e}")
            # Fallback to simple structural features
            features = np.array([
                structure.density,
                structure.volume / structure.num_sites,
                np.mean(structure.lattice.abc),
                np.std(structure.lattice.abc),
                structure.lattice.alpha,
                structure.lattice.beta,
                structure.lattice.gamma
            ])
        
        return features
    
    def fit_transform(self, features: np.ndarray) -> np.ndarray:
        """Fit scaler and transform features."""
        scaled_features = self.scaler.fit_transform(features)
        self.fitted = True
        return scaled_features
    
    def transform(self, features: np.ndarray) -> np.ndarray:
        """Transform features using fitted scaler."""
        if not self.fitted:
            raise ValueError("Must fit scaler first")
        return self.scaler.transform(features)


class RBTGraphNet(nn.Module):
    """Graph Neural Network for superconductor Tc prediction using RBT features."""
    
    def __init__(self, 
                 input_dim: int,
                 hidden_dim: int = 128,
                 num_layers: int = 3,
                 dropout: float = 0.2):
        super().__init__()
        
        if not TORCH_GEOMETRIC_AVAILABLE:
            raise ImportError("torch_geometric required for RBTGraphNet")
        
        self.num_layers = num_layers
        self.dropout = dropout
        
        # Graph convolution layers
        self.convs = nn.ModuleList()
        self.convs.append(GCNConv(input_dim, hidden_dim))
        
        for _ in range(num_layers - 1):
            self.convs.append(GCNConv(hidden_dim, hidden_dim))
        
        # Attention layer for important features
        self.attention = GATConv(hidden_dim, hidden_dim // 4, heads=4, concat=True)
        
        # Final prediction layers
        self.predictor = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, hidden_dim // 4),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 4, 1)
        )
    
    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        
        # Graph convolutions with residual connections
        for i, conv in enumerate(self.convs):
            if i == 0:
                x = F.relu(conv(x, edge_index))
            else:
                x_new = F.relu(conv(x, edge_index))
                x = x + x_new  # Residual connection
            x = F.dropout(x, p=self.dropout, training=self.training)
        
        # Attention mechanism
        x = self.attention(x, edge_index)
        x = F.dropout(x, p=self.dropout, training=self.training)
        
        # Global pooling
        x_mean = global_mean_pool(x, batch)
        x_max = global_max_pool(x, batch)
        x = torch.cat([x_mean, x_max], dim=1)
        
        # Final prediction
        out = self.predictor(x)
        return out.squeeze()


class RBTMLPredictor:
    """Complete ML prediction pipeline for superconductor discovery."""
    
    def __init__(self, model_type: str = 'random_forest'):
        self.model_type = model_type
        self.model = None
        self.feature_extractor = RBTFeatureExtractor()
        self.label_encoder = LabelEncoder()
        self.trained = False
        
        if model_type == 'random_forest':
            self.model = RandomForestRegressor(
                n_estimators=200,
                max_depth=20,
                min_samples_split=5,
                min_samples_leaf=2,
                random_state=42,
                n_jobs=-1
            )
        elif model_type == 'graph_net' and TORCH_GEOMETRIC_AVAILABLE:
            # Will be initialized during training
            pass
        else:
            raise ValueError(f"Unsupported model type: {model_type}")
    
    def prepare_training_data(self, 
                            supercon_data: List[Dict],
                            rbt_features: List[Dict]) -> Tuple[np.ndarray, np.ndarray]:
        """Prepare training data from SuperCon database and RBT features."""
        
        # Match compositions between datasets
        supercon_formulas = {entry['formula']: entry['tc_k'] for entry in supercon_data 
                           if 'tc_k' in entry and entry['tc_k'] > 0}
        
        X, y = [], []
        
        for rbt_entry in rbt_features:
            formula = rbt_entry.get('reduced_formula', '')
            
            # Look for matching superconductor data
            tc_value = None
            for sc_formula, tc in supercon_formulas.items():
                try:
                    if Composition(formula).reduced_formula == Composition(sc_formula).reduced_formula:
                        tc_value = tc
                        break
                except:
                    continue
            
            if tc_value is not None and tc_value > 0:
                # Extract features
                feature_vector = self._extract_feature_vector(rbt_entry)
                if feature_vector is not None:
                    X.append(feature_vector)
                    y.append(tc_value)
        
        print(f"Prepared {len(X)} training samples")
        return np.array(X), np.array(y)
    
    def _extract_feature_vector(self, rbt_entry: Dict) -> Optional[np.ndarray]:
        """Extract feature vector from RBT entry."""
        try:
            # RBT scores
            rbt_features = [
                rbt_entry.get('best_rbt_score', 0),
                rbt_entry.get('best_tau_valence', 0),
                rbt_entry.get('best_kappa_curvature', 0),
                rbt_entry.get('best_phase_locking', 0),
                rbt_entry.get('structure_quality', 0.5)
            ]
            
            # Composition features
            comp_features = [
                rbt_entry.get('n_elements', 0),
                rbt_entry.get('n_atoms', 0),
                rbt_entry.get('metal_fraction', 0),
                rbt_entry.get('avg_electronegativity', 0),
                rbt_entry.get('electronegativity_spread', 0),
                rbt_entry.get('avg_mass', 0),
                rbt_entry.get('density_estimate', 0)
            ]
            
            # Combine all features
            features = np.array(rbt_features + comp_features)
            
            # Check for invalid values
            if np.any(np.isnan(features)) or np.any(np.isinf(features)):
                return None
                
            return features
            
        except Exception as e:
            print(f"Warning: Feature extraction failed for {rbt_entry.get('formula', 'unknown')}: {e}")
            return None
    
    def train(self, X: np.ndarray, y: np.ndarray, validation_split: float = 0.2) -> Dict[str, float]:
        """Train the ML model."""
        
        # Split data
        X_train, X_val, y_train, y_val = train_test_split(
            X, y, test_size=validation_split, random_state=42
        )
        
        # Scale features
        X_train_scaled = self.feature_extractor.fit_transform(X_train)
        X_val_scaled = self.feature_extractor.transform(X_val)
        
        # Train model
        if self.model_type == 'random_forest':
            self.model.fit(X_train_scaled, y_train)
            
            # Predictions
            y_train_pred = self.model.predict(X_train_scaled)
            y_val_pred = self.model.predict(X_val_scaled)
            
        elif self.model_type == 'graph_net':
            # Would implement graph net training here
            raise NotImplementedError("Graph net training not implemented yet")
        
        # Calculate metrics
        train_metrics = {
            'train_mae': mean_absolute_error(y_train, y_train_pred),
            'train_rmse': np.sqrt(mean_squared_error(y_train, y_train_pred)),
            'train_r2': r2_score(y_train, y_train_pred)
        }
        
        val_metrics = {
            'val_mae': mean_absolute_error(y_val, y_val_pred),
            'val_rmse': np.sqrt(mean_squared_error(y_val, y_val_pred)),
            'val_r2': r2_score(y_val, y_val_pred)
        }
        
        self.trained = True
        
        metrics = {**train_metrics, **val_metrics}
        return metrics
    
    def predict(self, X: np.ndarray) -> np.ndarray:
        """Make predictions on new data."""
        if not self.trained:
            raise ValueError("Model must be trained first")
        
        X_scaled = self.feature_extractor.transform(X)
        return self.model.predict(X_scaled)
    
    def predict_candidates(self, candidates: List[Dict]) -> List[Dict]:
        """Predict Tc for candidate materials."""
        if not self.trained:
            raise ValueError("Model must be trained first")
        
        # Extract features
        X = []
        valid_candidates = []
        
        for candidate in candidates:
            feature_vector = self._extract_feature_vector(candidate)
            if feature_vector is not None:
                X.append(feature_vector)
                valid_candidates.append(candidate)
        
        if not X:
            print("Warning: No valid candidates for prediction")
            return []
        
        # Make predictions
        X = np.array(X)
        predictions = self.predict(X)
        
        # Add predictions to candidates
        for candidate, tc_pred in zip(valid_candidates, predictions):
            candidate['ml_tc_prediction'] = float(tc_pred)
            candidate['ml_confidence'] = self._estimate_confidence(candidate)
        
        return valid_candidates
    
    def _estimate_confidence(self, candidate: Dict) -> float:
        """Estimate prediction confidence (placeholder)."""
        # This would use uncertainty quantification in a real implementation
        rbt_score = candidate.get('best_rbt_score', 0)
        quality = candidate.get('structure_quality', 0.5)
        return min(rbt_score * quality, 1.0)
    
    def save_model(self, filepath: str):
        """Save trained model."""
        model_data = {
            'model': self.model,
            'feature_extractor': self.feature_extractor,
            'model_type': self.model_type,
            'trained': self.trained
        }
        
        with open(filepath, 'wb') as f:
            pickle.dump(model_data, f)
    
    def load_model(self, filepath: str):
        """Load trained model."""
        with open(filepath, 'rb') as f:
            model_data = pickle.load(f)
        
        self.model = model_data['model']
        self.feature_extractor = model_data['feature_extractor']
        self.model_type = model_data['model_type']
        self.trained = model_data['trained']


def main():
    parser = argparse.ArgumentParser(
        description="ML prediction of superconducting Tc using RBT features"
    )
    
    # Mode selection
    parser.add_argument(
        "--train", action="store_true",
        help="Train mode: build ML model from data"
    )
    parser.add_argument(
        "--predict", action="store_true", 
        help="Predict mode: use trained model on new candidates"
    )
    parser.add_argument(
        "--rank", type=str,
        help="Rank candidates by ML predictions (provide candidates file)"
    )
    
    # Training arguments
    parser.add_argument(
        "--supercon_data", type=str,
        help="Path to SuperCon database CSV file"
    )
    parser.add_argument(
        "--rbt_features", type=str,
        help="Path to RBT features CSV file"
    )
    parser.add_argument(
        "--model_type", choices=['random_forest', 'graph_net'], default='random_forest',
        help="Type of ML model to use"
    )
    
    # Prediction arguments
    parser.add_argument(
        "--model", type=str,
        help="Path to trained model file"
    )
    parser.add_argument(
        "--candidates", type=str,
        help="Path to candidates file for prediction"
    )
    parser.add_argument(
        "--top", type=int, default=200,
        help="Number of top candidates to output"
    )
    
    # Output arguments
    parser.add_argument(
        "--out", type=str, default="results/ml_predictions.csv",
        help="Output file path"
    )
    parser.add_argument(
        "--model_out", type=str, default="models/best_model.pkl",
        help="Path to save trained model"
    )
    
    args = parser.parse_args()
    
    # Create output directories
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    if args.model_out:
        Path(args.model_out).parent.mkdir(parents=True, exist_ok=True)
    
    # Initialize predictor
    predictor = RBTMLPredictor(model_type=args.model_type)
    
    if args.train:
        # Training mode
        if not args.supercon_data or not args.rbt_features:
            raise ValueError("Training requires --supercon_data and --rbt_features")
        
        print("Loading training data...")
        
        # Load SuperCon data
        supercon_data = load_supercon_data(args.supercon_data)
        print(f"Loaded {len(supercon_data)} SuperCon entries")
        
        # Load RBT features  
        rbt_df = pd.read_csv(args.rbt_features)
        rbt_features = rbt_df.to_dict('records')
        print(f"Loaded {len(rbt_features)} RBT feature entries")
        
        # Prepare training data
        X, y = predictor.prepare_training_data(supercon_data, rbt_features)
        
        if len(X) < 10:
            print("Warning: Very few training samples. Results may be unreliable.")
        
        # Train model
        print(f"Training {args.model_type} model on {len(X)} samples...")
        metrics = predictor.train(X, y)
        
        # Print results
        print("\nTraining Results:")
        for metric, value in metrics.items():
            print(f"{metric}: {value:.4f}")
        
        # Save model
        predictor.save_model(args.model_out)
        print(f"Saved model to {args.model_out}")
        
        # Feature importance (for random forest)
        if args.model_type == 'random_forest':
            feature_names = [
                'rbt_score', 'tau_valence', 'kappa_curvature', 'phase_locking', 'structure_quality',
                'n_elements', 'n_atoms', 'metal_fraction', 'avg_electronegativity', 
                'electronegativity_spread', 'avg_mass', 'density_estimate'
            ]
            
            importance = predictor.model.feature_importances_
            
            print("\nFeature Importance:")
            for name, imp in zip(feature_names, importance):
                print(f"{name}: {imp:.4f}")
    
    elif args.predict or args.rank:
        # Prediction mode
        if not args.model:
            raise ValueError("Prediction requires --model")
        
        print(f"Loading model from {args.model}...")
        predictor.load_model(args.model)
        
        # Load candidates
        candidates_file = args.candidates or args.rank
        if candidates_file.endswith('.csv'):
            candidates_df = pd.read_csv(candidates_file)
            candidates = candidates_df.to_dict('records')
        else:
            with open(candidates_file, 'r') as f:
                candidates = json.load(f)
        
        print(f"Making predictions for {len(candidates)} candidates...")
        
        # Make predictions
        predicted_candidates = predictor.predict_candidates(candidates)
        
        # Sort by ML prediction
        predicted_candidates.sort(key=lambda x: x.get('ml_tc_prediction', 0), reverse=True)
        
        # Take top candidates
        top_candidates = predicted_candidates[:args.top]
        
        # Save results
        df_out = pd.DataFrame(top_candidates)
        
        # Select key columns
        key_columns = [
            'reduced_formula', 'ml_tc_prediction', 'ml_confidence',
            'best_rbt_score', 'best_tc_estimate', 'best_tau_valence',
            'best_kappa_curvature', 'best_phase_locking'
        ]
        
        available_columns = [col for col in key_columns if col in df_out.columns]
        df_output = df_out[available_columns]
        
        df_output.to_csv(args.out, index=False, float_format='%.4f')
        
        print(f"\nSaved {len(top_candidates)} predictions to {args.out}")
        
        # Print top 10 predictions
        print("\nTop 10 ML predictions:")
        print("Rank | Formula | ML Tc | Confidence | RBT Score")
        print("-" * 60)
        
        for i, candidate in enumerate(top_candidates[:10]):
            print(f"{i+1:4d} | {candidate.get('reduced_formula', 'N/A'):12s} | "
                  f"{candidate.get('ml_tc_prediction', 0):6.1f}K | "
                  f"{candidate.get('ml_confidence', 0):8.3f} | "
                  f"{candidate.get('best_rbt_score', 0):8.4f}")
    
    else:
        print("Must specify either --train or --predict/--rank mode")


if __name__ == "__main__":
    main() 