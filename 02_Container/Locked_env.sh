echo ""
echo "Exporting locked conda environments..."


CONTAINER_DIR="02_Container"


# No conda environments created
if [ ! -d ".snakemake/conda" ] || [ -z "$(ls -A .snakemake/conda)" ]; then
    echo "No conda environments were created."
    echo "Nothing to export."
    exit 0
fi


# Get original yaml files only
yaml_files=(
    $(find "$CONTAINER_DIR" \
        -maxdepth 1 \
        -name "*.yaml" \
        ! -name "*_locked.yaml" \
        | sort)
)


# Loop over Snakemake environments
for env in .snakemake/conda/*_; do

    if [ ! -d "$env" ]; then
        continue
    fi


    env_hash=$(basename "$env")

    echo ""
    echo "Checking environment: $env_hash"


    matched=false


    # Snakemake stores the yaml used for this environment here
    snakemake_yaml=".snakemake/conda/${env_hash}.yaml"


    if [ ! -f "$snakemake_yaml" ]; then
        echo "No yaml found for $env_hash"
        continue
    fi


    # Compare against original yaml files
    for yaml in "${yaml_files[@]}"; do

        name=$(basename "$yaml" .yaml)


        if diff -q "$yaml" "$snakemake_yaml" >/dev/null; then

            echo "Matched: $name"


            conda env export \
                --prefix "$env" \
                --no-builds \
                > "${CONTAINER_DIR}/${name}_locked.yaml"


            echo "Created:"
            echo "  ${CONTAINER_DIR}/${name}_locked.yaml"


            matched=true
            break

        fi

    done


    # Fallback if no match was found
    if [ "$matched" = false ]; then

        echo "WARNING: Could not identify $env_hash"
        echo "Exporting with generic name."


        conda env export \
            --prefix "$env" \
            --no-builds \
            > "${CONTAINER_DIR}/${env_hash}_locked.yaml"


        echo "Created:"
        echo "  ${CONTAINER_DIR}/${env_hash}_locked.yaml"

    fi

done


echo ""
echo "Locked environments exported successfully!"