# uuidgen
rm -rf transcriptomics.55fe0636-c9f9-4300-852e-14b3b9c86f3b.zip
zip -r transcriptomics.55fe0636-c9f9-4300-852e-14b3b9c86f3b.zip test.sh testData/ nextflow.config  main.nf nextflow  installDocker.sh README.md

ossutil64 cp -r transcriptomics.55fe0636-c9f9-4300-852e-14b3b9c86f3b.zip  oss://wangyang-bucket/pay/ --force