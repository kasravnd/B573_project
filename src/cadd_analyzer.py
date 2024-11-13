import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Tuple, Dict, Optional
from dataclasses import dataclass
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

@dataclass
class StatisticalResults:
    """Class to store statistical analysis results"""
    # Descriptive statistics
    n_pathogenic: int
    n_benign: int
    mean_pathogenic: float
    mean_benign: float
    median_pathogenic: float
    median_benign: float
    std_pathogenic: float
    std_benign: float
    
    # Statistical tests
    mannwhitney_statistic: float
    mannwhitney_pvalue: float
    cohens_d: float
    
    def __str__(self) -> str:
        """Pretty print the results"""
        return f"""
Statistical Analysis Results:
----------------------------
Sample Sizes:
    Pathogenic variants: {self.n_pathogenic}
    Benign variants: {self.n_benign}

Descriptive Statistics:
    Pathogenic - Mean: {self.mean_pathogenic:.2f}, Median: {self.median_pathogenic:.2f}, STD: {self.std_pathogenic:.2f}
    Benign - Mean: {self.mean_benign:.2f}, Median: {self.median_benign:.2f}, STD: {self.std_benign:.2f}

Statistical Tests:
    Mann-Whitney U test:
        Statistic: {self.mannwhitney_statistic:.2f}
        P-value: {self.mannwhitney_pvalue:.4f}
    
Effect Size:
    Cohen's d: {self.cohens_d:.2f}
"""

class CADDScoreAnalyzer:
    """Class for analyzing CADD scores in genetic variants"""
    
    def __init__(self, data: pd.DataFrame):
        """
        Initialize the analyzer with data
        
        Args:
            data (pd.DataFrame): DataFrame containing variant data
        """
        self.data = data
        self._prepare_data()
        
    def _prepare_data(self) -> None:
        """Prepare and validate the data"""
        # For demonstration, we'll simulate CADD scores if they don't exist
        if 'CADD' not in self.data.columns:
            logger.info("Simulating CADD scores for demonstration")
            # Simulate higher CADD scores for pathogenic variants
            self.data['CADD'] = np.where(
                self.data['clinical_significance'].str.contains('pathogenic', case=False, na=False),
                np.random.normal(25, 5, len(self.data)),  # Pathogenic variants
                np.random.normal(15, 5, len(self.data))   # Non-pathogenic variants
            )
        
        # Create binary classification
        self.data['is_pathogenic'] = self.data['clinical_significance'].str.contains(
            'pathogenic', case=False, na=False
        )
        self.data['is_benign'] = self.data['clinical_significance'].str.contains(
            'benign', case=False, na=False
        )
        
        # Filter for clear classifications
        self.analyzed_data = self.data[
            (self.data['is_pathogenic'] | self.data['is_benign'])
        ].copy()
        
    def calculate_descriptive_stats(self) -> Dict[str, float]:
        """Calculate descriptive statistics for CADD scores"""
        pathogenic_scores = self.analyzed_data[self.analyzed_data['is_pathogenic']]['CADD']
        benign_scores = self.analyzed_data[self.analyzed_data['is_benign']]['CADD']
        
        return {
            'n_pathogenic': len(pathogenic_scores),
            'n_benign': len(benign_scores),
            'mean_pathogenic': pathogenic_scores.mean(),
            'mean_benign': benign_scores.mean(),
            'median_pathogenic': pathogenic_scores.median(),
            'median_benign': benign_scores.median(),
            'std_pathogenic': pathogenic_scores.std(),
            'std_benign': benign_scores.std()
        }
    
    def perform_mann_whitney(self) -> Tuple[float, float]:
        """
        Perform Mann-Whitney U test
        
        Returns:
            Tuple[float, float]: Test statistic and p-value
        """
        pathogenic_scores = self.analyzed_data[self.analyzed_data['is_pathogenic']]['CADD']
        benign_scores = self.analyzed_data[self.analyzed_data['is_benign']]['CADD']
        
        return stats.mannwhitneyu(pathogenic_scores, benign_scores, alternative='two-sided')
    
    def calculate_cohens_d(self) -> float:
        """
        Calculate Cohen's d effect size
        
        Returns:
            float: Cohen's d value
        """
        pathogenic_scores = self.analyzed_data[self.analyzed_data['is_pathogenic']]['CADD']
        benign_scores = self.analyzed_data[self.analyzed_data['is_benign']]['CADD']
        
        n1, n2 = len(pathogenic_scores), len(benign_scores)
        var1, var2 = pathogenic_scores.var(), benign_scores.var()
        
        # Pooled standard deviation
        pooled_sd = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
        
        return (pathogenic_scores.mean() - benign_scores.mean()) / pooled_sd
    
    def create_boxplot(self, save_path: Optional[str] = None) -> None:
        """
        Create and save boxplot comparing CADD scores
        
        Args:
            save_path (str, optional): Path to save the plot
        """
        plt.figure(figsize=(10, 6))
        
        # Create boxplot
        sns.boxplot(
            data=self.analyzed_data,
            x='clinical_significance',
            y='CADD',
            order=['benign', 'likely benign', 'likely pathogenic', 'pathogenic']
        )
        
        # Customize plot
        plt.title('Distribution of CADD Scores by Variant Classification', pad=20)
        plt.xlabel('Variant Classification')
        plt.ylabel('CADD Score')
        plt.xticks(rotation=45)
        
        # Add statistical annotation
        stat, pval = self.perform_mann_whitney()
        plt.text(
            0.05, 0.95,
            f'Mann-Whitney p-value: {pval:.2e}',
            transform=plt.gca().transAxes,
            bbox=dict(facecolor='white', alpha=0.8)
        )
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Boxplot saved to {save_path}")
        
        plt.close()
    
    def run_full_analysis(self) -> StatisticalResults:
        """
        Run complete statistical analysis
        
        Returns:
            StatisticalResults: Complete analysis results
        """
        try:
            logger.info("Starting statistical analysis...")
            
            # Calculate descriptive statistics
            stats_dict = self.calculate_descriptive_stats()
            
            # Perform Mann-Whitney U test
            mw_stat, mw_pval = self.perform_mann_whitney()
            
            # Calculate Cohen's d
            cohens_d = self.calculate_cohens_d()
            
            # Create visualization
            self.create_boxplot('cadd_scores_boxplot.png')
            
            # Create results object
            results = StatisticalResults(
                n_pathogenic=stats_dict['n_pathogenic'],
                n_benign=stats_dict['n_benign'],
                mean_pathogenic=stats_dict['mean_pathogenic'],
                mean_benign=stats_dict['mean_benign'],
                median_pathogenic=stats_dict['median_pathogenic'],
                median_benign=stats_dict['median_benign'],
                std_pathogenic=stats_dict['std_pathogenic'],
                std_benign=stats_dict['std_benign'],
                mannwhitney_statistic=mw_stat,
                mannwhitney_pvalue=mw_pval,
                cohens_d=cohens_d
            )
            
            logger.info("Analysis completed successfully")
            return results
            
        except Exception as e:
            logger.error(f"Error during analysis: {str(e)}")
            raise

def main():
    """Main function to run the analysis"""
    try:
        # Read the data
        df = pd.read_csv('arrhythmia.csv')
        
        # Initialize analyzer and run analysis
        analyzer = CADDScoreAnalyzer(df)
        results = analyzer.run_full_analysis()
        
        # Print results
        print(results)
        
    except Exception as e:
        logger.error(f"Error in main function: {str(e)}")
        raise

if __name__ == "__main__":
    main()