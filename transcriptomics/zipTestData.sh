# uuidgen
rm -rf transcriptomics.test.data.55fe0636-c9f9-4300-852e-14b3b9c86f3b.zip
zip -r transcriptomics.test.data.55fe0636-c9f9-4300-852e-14b3b9c86f3b.zip testData/

ossutil64 cp -r transcriptomics.test.data.55fe0636-c9f9-4300-852e-14b3b9c86f3b.zip  oss://wangyang-bucket/pay/ --force