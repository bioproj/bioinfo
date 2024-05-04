# .libPaths(c("/XXXX/R_LIBS",.libPaths()))
.libPaths()
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) 
options()$repos  
options()$BioC_mirror 

BiocManager::version()