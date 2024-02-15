import pandas as pd

from snakemake.io import Wildcards

class Samples:
    """
    Convert the OTP-exported metadata spreadsheet into a pandas.DataFrame
    that provides metadata for the workflow.
    """
               
    def __init__(self, config):
        IDENTIFIERS = config["metadata"]["identifiers"]
        metadata_files = config["metadata"]["table"]

        self.output_base_dir = config["paths"]["output_dir"]
        self.target_templates = config["paths"]["target_templates"]
        
        # columns to select for subsampling 
        columns = ["PID",
                   "Sample Type",
                   "Species Scientific Name",
                   "STRAIN",
                   "FastQ Path",
                   "FastQC Path",
                   "Run ID",
                   "ILSE_NO",
                   "TISSUE",
                   "BIRTH",
                   "DATE_OF_BIRTH",
                   "DATE_OF_DEATH",
                   "LANE_NO",
                   "READ",
                   "CELLRANGER_FASTQ_PATH",
                   "individual"]  # include 'individual' col to be generated below

        metadata_full = pd.concat((pd.read_csv(f) for f in metadata_files), ignore_index=True)

        # if required, generate 'individual' column containing all concatenated identifiers
        if not "individual" in metadata_full.columns: 
            metadata_full['individual'] = metadata_full[IDENTIFIERS].apply('_'.join, axis=1)
   
        metadata_full = self.get_cellranger_filename(metadata_full)

        self.metadata = self.select_columns(metadata_full, custom_columns = columns + IDENTIFIERS)

    def get_cellranger_filename(self, df: pd.DataFrame) -> pd.DataFrame:

        """Add column containing CellRanger compatible filename,
        i.e. in the format of
        [Sample Name]_S[Sample_Number]_L00[Lane Number]_[Read Type]_001.fastq.gz
        
        Here, [Sample Name] consists of "individual".
        
        See also: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input
        """

        grouped = df.groupby(["individual", "LANE_NO", "READ"])
        groups = []
        for name, group in grouped:
            group = group.sort_values("FastQ Path")  # import to have consistent sorting
            group["multi_sample_idx"] = range(1, len(group)+1)
            group["CELLRANGER_FASTQ_PATH"] = group.agg(
                "{0[individual]}_S{0[multi_sample_idx]}_L00{0[LANE_NO]}_R{0[READ]}_001.fastq.gz".format,
                axis=1,
            )
            groups.append(group)

        return pd.concat(groups)

    def select_columns(self,
                       df: pd.DataFrame,
                       custom_columns: list = None):
        """Select/Subset columns from DataFrame to reduce
        DataFrame dimensions. A list of column names can be provided by `custom_columns` """
        
        if custom_columns:
          print(custom_columns)
          df_subset = df[custom_columns] 
        else:
          print(custom_columns)
          df_subset = df[columns]
        
        return df_subset
            
    @staticmethod
    def filter_by_wildcards(
            wildcards: Wildcards,
            data: pd.DataFrame,
            column: str = None,
            unique: bool = True,
            **add_filters,
    ) -> list:
        """Get item from DataFrame by subsetting using the following index
        attributes provided via Snakemakes wildcards object.
        :return param:
        """

        _identifiers = ["individual", "Species_ID"] 

        if wildcards:
            filters = dict((k, getattr(wildcards, k)) for k in _identifiers)
        else:
            filters = {}

        filters.update(add_filters)
        subset_rows = (
            (data.reset_index()[list(filters)] == pd.Series(filters)).all(axis=1).to_list()
        )  # noqa: E501

        params = data.loc[subset_rows].get(column, None)  # select column `fieldname`
        if params is None:
            return [None]
        elif unique:
            return params.drop_duplicates().to_list()
        else:
            return params.to_list()

    @property
    def targets(self) -> list:
        """Create all target for the workflow."""

        output_base_path = self.output_base_dir

        target_templates = dict((cat, f"{output_base_path}/{path}") for cat, path in self.target_templates.items())

        targets = list(map(lambda x: self.metadata.agg(x.format, axis=1).drop_duplicates().to_list(), target_templates.values()))
        return targets
      
      
      
#grouped = metadata_subset.groupby(["PID", "Sample Type"])["FastQ Path"].apply(list)
