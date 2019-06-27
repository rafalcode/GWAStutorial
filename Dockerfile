FROM rocker/r-ver:3.5.3
# The rocker project has a number of Debian based docker images for
# several versions of R. When looking to upgrade, look at
# https://github.com/rocker-org/rocker

# ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive
RUN echo "deb http://deb.debian.org/debian stretch-backports main" >> /etc/apt/sources.list \
    && sed -i 's/main/main contrib non-free/g' /etc/apt/sources.list \
    && apt-get update \
    && apt-get install -y apt-utils \
    && apt-get install -y build-essential \
    && apt-get install -y python \
    && apt-get install -y python3 \
    && apt-get install -y python3-pandas python3-pandas-lib python3-numpy \
    && apt-get install -y wget git \
    && apt-get install -y bzip2 tree \
    && apt-get install -y zlib1g-dev libbz2-1.0 libbz2-dev liblzma5 liblzma-dev \
    && apt-get install -y vim pkg-config libnlopt-dev \
    && apt-get install -y libssl-dev libcurl4-openssl-dev plink1.9

RUN Rscript -e "update.packages(ask=FALSE)" \
    && Rscript -e "install.packages('logging')" \
    && Rscript -e "install.packages('vcfR')" \
    && Rscript -e "install.packages('doParallel')" \
    && Rscript -e "install.packages('plyr')" \
    && Rscript -e "install.packages('dplyr')" \
    && Rscript -e "install.packages('randomForest')" \
    && Rscript -e "install.packages('xgboost')" \
    && Rscript -e "install.packages('stringr')" \
    && Rscript -e "install.packages('glmnet')" \
    && Rscript -e "install.packages('devtools')" \
    && Rscript -e "install.packages('car')" \
    && Rscript -e "library(devtools); install_github('rafalcode/GenABEL.data'); install_github('rafalcode/GenABEL');" \
    && Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite(); biocLite('snpStats'); biocLite('SNPRelate')"

#hand compile bftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 \
    && tar xjf bcftools-1.9.tar.bz2 \
    && cd bcftools-1.9 \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && rm -rf bcftools-1.9.tar.bz2 bcftools-1.9

# and samtools, why not?
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && tar xjf samtools-1.9.tar.bz2 \
    && cd samtools-1.9/htslib-1.9 \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && ./configure --without-curses \
    && make \
    && make install \
    && cd .. \
    && rm -rf samtools-1.9 samtools-1.9.tar.bz2 
