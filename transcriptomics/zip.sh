# uuidgen
rm -rf transcriptomics.pay.3fcb3ca7-9d87-4195-ac7b-1eb7fdd75edc.zip
zip -r transcriptomics.pay.3fcb3ca7-9d87-4195-ac7b-1eb7fdd75edc.zip bin/ rstudio.sh env.R test.sh nextflow.config  main.nf nextflow  dockerPull.sh dockerInstall.sh README.md downlaodTestData.sh

ossutil64 cp -r transcriptomics.pay.3fcb3ca7-9d87-4195-ac7b-1eb7fdd75edc.zip  oss://wangyang-bucket/pay/ --force