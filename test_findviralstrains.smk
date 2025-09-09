# test_findviralstrains.smk
include: "findviralstrains.smk"

# Override or add test-specific rules
rule test_validation:
    input:
        # List all expected final outputs
        expand("output/{analysis_ID}/decomp_results.csv", analysis_ID=config["analysis_ID"]),
        expand("output/{analysis_ID}/visualizations/", analysis_ID=config["analysis_ID"])
    output:
        touch("test_success.txt")
    shell:
        """
        # Add validation commands here
        echo "Pipeline completed successfully, validating outputs..."
        
        # Example validations
        python scripts/validate_outputs.py {input[0]}
        """
