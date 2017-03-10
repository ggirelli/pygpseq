#
#
# ./pyGPSeq_run_batch.sh paramTable
#
# Run all the pyGPSeq jobs specified in the provided paramTable.
# 


# CHECK PARAMETERS
# ========================

# Check number of parameters
if [ "$#" -ne 1 ]; then
	echo -e "Correct usage: ./pyGPSeq_run_batch.sh paramTable";
	exit 1;
fi

# Check if the paramTable exists
paramTable=$1
if [ ! -f $paramTable ]; then
	echo -e "Correct usage: ./pyGPSeq_run_batch.sh paramTable";
	echo -e "paramTable not found at '$paramTable'.\nSTOPPED";
	exit 1;
fi

# RUN
# =======================

# Count rows in paramTable
nrow=`wc -l $paramTable | cut -d ' ' -f 1`;

# Remove 1 for the header
nrow=$(($nrow - 1));

# Cycle through the rows
for i in $(seq 1 $nrow); do
	./pyGPSeq_run_from_table.py $paramTable $i
done

