LOGGING_FMT = "%(name)-20s %(levelname)-7s @ %(asctime)s: %(message)s"
LOGGING_DATE_FMT = "%m/%d/%y %H:%M:%S"

DEFAULT_PIPELINE_CONFIG = {
    "global": {
        "number_cpus": 8,
        "working_dir": "./",
        "log_level": 20,
    },
    "data": {
        "input_dir": None,
        "fasta": None,
        "bwa_index_prefix": None,
        "restriction_site": "T^TAA",
        "restriction_name": "MseI",
        "chromosome_file": "hg19",
    },
    "processes": {
        "keep": [],
        "is_qc": True,
    },
    "extract_pet": {
        "adapter": "auto",
        "mismatch_adapter": 0.3,
        "n_fq_for_infer": 1000,
        "search_start_pos": 70,
        "prob_thresh": 0.75,
        "linker_A": "",
        "linker_B": "",
        "mismatch": 0.3,
        "pet_len_range": (10, 22),
        "pet_cut_len": 20,
    },
    "build_bedpe": {
        "mapq": 20,
        "iterations": 7,
    },
    "noise_reduce": {
        "restriction_sites_file": None,
        "threshold_span": -1,
    },
    "quality_control": {
        "report_format": "html",
        "long_range_cutoff": 5000,
    },
    "result": {
        "result_formats": ['.cool'],
        "resolutions": [2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000],
        "juicer_tools_jar": None,
    },
}
