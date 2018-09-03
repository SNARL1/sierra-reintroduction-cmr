# Load dependencies
library(rstan)
library(lubridate)
library(ggthemes)
library(reshape2)
library(assertthat)
library(tidyverse)
library(mgcv)
library(ggridges)
library(patchwork)
library(scales)
library(ggrepel)
theme_set(theme_minimal() + 
            theme(panel.grid.minor = element_blank(), 
                  panel.grid.major = element_line(colour = 'grey90')))
source('R/helpers.R')
