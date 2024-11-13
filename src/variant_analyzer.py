import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple
import logging
from collections import Counter

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ArrhythmiaAnalyzer:
    """Class for analyzing genetic variants in cardiac arrhythmia"""
    
    def __init__(self, data: pd.DataFrame):
        """Initialize with DataFrame"""
        self.data = data
        self._clean_data()
        
    def _clean_data(self) -> None:
        """Clean and prepare the data"""
        # Remove rows where clinical_significance is NA
        self.data = self.data.dropna(subset=['clinical_significance'])
        
        # Standardize clinical significance categories
        self.data['clinical_significance'] = self.data['clinical_significance'].str.lower()
        
    def analyze_variant_classifications(self) -> Dict[str, int]:
        """Analyze distribution of variant classifications"""
        return dict(self.data['clinical_significance'].value_counts())
    
    def analyze_genes(self) -> Dict[str, int]:
        """Analyze frequency of variants in different genes"""
        return dict(self.data['attributes_associated_gene'].value_counts().head(10))
    
    def analyze_variant_types(self) -> Dict[str, int]:
        """Analyze distribution of variant types"""
        return dict(self.data['var_class'].value_counts())
    
    def analyze_population_frequencies(self) -> pd.DataFrame:
        """Analyze variant frequencies across different populations"""
        pop_columns = [col for col in self.data.columns if col.startswith('pop.gnomad_')]
        return self.data[pop_columns].agg(['mean', 'median', 'std']).round(4)
    
    def create_variant_classification_plot(self, save_path: str = None) -> None:
        """Create bar plot of variant classifications"""
        plt.figure(figsize=(12, 6))
        sns.countplot(data=self.data, y='clinical_significance', order=self.data['clinical_significance'].value_counts().index)
        plt.title('Distribution of Variant Classifications')
        plt.xlabel('Count')
        plt.ylabel('Clinical Significance')
        
        if save_path:
            plt.savefig(save_path, bbox_inches='tight')
        plt.close()
    
    def create_gene_distribution_plot(self, save_path: str = None) -> None:
        """Create bar plot of variant distribution across genes"""
        plt.figure(figsize=(12, 6))
        gene_counts = self.data['attributes_associated_gene'].value_counts().head(10)
        sns.barplot(x=gene_counts.values, y=gene_counts.index)
        plt.title('Top 10 Genes with Most Variants')
        plt.xlabel('Number of Variants')
        plt.ylabel('Gene')
        
        if save_path:
            plt.savefig(save_path, bbox_inches='tight')
        plt.close()
    
    def create_population_frequency_heatmap(self, save_path: str = None) -> None:
        """Create heatmap of variant frequencies across populations"""
        pop_columns = [col for col in self.data.columns if col.startswith('pop.gnomad_')]
        plt.figure(figsize=(12, 8))
        
        # Calculate correlation matrix
        corr_matrix = self.data[pop_columns].corr()
        
        # Create heatmap
        sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0)
        plt.title('Correlation of Variant Frequencies Across Populations')
        
        if save_path:
            plt.savefig(save_path, bbox_inches='tight')
        plt.close()
    
    def analyze_pathogenic_variants(self) -> pd.DataFrame:
        """Analyze characteristics of pathogenic variants"""
        pathogenic_data = self.data[self.data['clinical_significance'].str.contains('pathogenic')]
        
        analysis = {
            'total_count': len(pathogenic_data),
            'genes_affected': pathogenic_data['attributes_associated_gene'].nunique(),
            'variant_types': dict(pathogenic_data['var_class'].value_counts()),
            'consequences': dict(pathogenic_data['most_severe_consequence'].value_counts().head(5))
        }
        
        return pd.DataFrame([analysis])
    
    def run_full_analysis(self) -> Dict:
        """Run complete analysis pipeline"""
        try:
            logger.info("Starting analysis...")
            
            results = {
                'variant_classifications': self.analyze_variant_classifications(),
                'gene_distribution': self.analyze_genes(),
                'variant_types': self.analyze_variant_types(),
                'population_frequencies': self.analyze_population_frequencies(),
                'pathogenic_analysis': self.analyze_pathogenic_variants()
            }
            
            # Create visualizations
            self.create_variant_classification_plot('variant_classifications.png')
            self.create_gene_distribution_plot('gene_distribution.png')
            self.create_population_frequency_heatmap('population_frequencies.png')
            
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
        analyzer = ArrhythmiaAnalyzer(df)
        results = analyzer.run_full_analysis()
        
        # Print results
        print("\nAnalysis Results:")
        print("\nVariant Classifications:")
        for classification, count in results['variant_classifications'].items():
            print(f"{classification}: {count}")
            
        print("\nTop 5 Genes with Most Variants:")
        for gene, count in list(results['gene_distribution'].items())[:5]:
            print(f"{gene}: {count}")
            
        print("\nVariant Types:")
        for var_type, count in results['variant_types'].items():
            print(f"{var_type}: {count}")
            
        print("\nPopulation Frequency Statistics:")
        print(results['population_frequencies'])
        
        print("\nPathogenic Variant Analysis:")
        print(results['pathogenic_analysis'])
        
    except Exception as e:
        logger.error(f"Error in main function: {str(e)}")
        raise

if __name__ == "__main__":
    main()