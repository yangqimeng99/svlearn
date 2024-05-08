inputFa=$1
outPrefix=$2

genmap index -F ${inputFa} -I ./`basename ${inputFa}`.index

genmap map \
	-K 50 -E 1 \
	-I ./`basename ${inputFa}`.index \
	-O ${outPrefix} \
	--threads 6 \
	-fl --txt
