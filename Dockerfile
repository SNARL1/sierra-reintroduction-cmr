FROM rocker/tidyverse:latest
MAINTAINER Max Joseph <maxwell.b.joseph@colorado.edu>

# rstan installation taken from jrnold/rstan
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
	apt-utils \
	ed \
	libnlopt-dev \
	pdftk \
	texlive-full \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/

# Global site-wide config -- neeeded for building packages
RUN mkdir -p $HOME/.R/ \
    && echo "CXXFLAGS=-O3 -mtune=native -march=native -Wno-unused-variable -Wno-unused-function -flto -ffat-lto-objects  -Wno-unused-local-typedefs \n" >> $HOME/.R/Makevars

# Config for rstudio user
RUN mkdir -p $HOME/.R/ \
    && echo "CXXFLAGS=-O3 -mtune=native -march=native -Wno-unused-variable -Wno-unused-function -flto -ffat-lto-objects  -Wno-unused-local-typedefs -Wno-ignored-attributes -Wno-deprecated-declarations\n" >> $HOME/.R/Makevars \
    && echo "rstan::rstan_options(auto_write = TRUE)\n" >> /home/rstudio/.Rprofile \
    && echo "options(mc.cores = parallel::detectCores())\n" >> /home/rstudio/.Rprofile

# Install rstan
RUN install2.r --error --deps TRUE rstan

# Install other packages
RUN install2.r --error --deps TRUE \
    assertthat \
    ggrepel \
    ggridges \
    ggthemes \
    knitcitations \
    lubridate \
    mgcv \
    reshape2 \
    scales

RUN installGithub.r thomasp85/patchwork \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

WORKDIR /home/rstudio/sierra-reintroduction-cmr

COPY . .

RUN chown rstudio . -R

