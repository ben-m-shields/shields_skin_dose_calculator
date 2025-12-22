FROM rocker/shiny-verse:latest

# ------------------------------------------------------------
# System dependencies
# ------------------------------------------------------------
RUN apt-get update && apt-get install -y \
    curl \
    git \
    wget \
    fonts-liberation \
    libcups2 \
    libgtk-3-0 \
    libxcomposite1 \
    libxfixes3 \
    xdg-utils \
    libnss3 \
    libatk-bridge2.0-0 \
    libxkbcommon0 \
    libxdamage1 \
    libxrandr2 \
    libgbm1 \
    libasound2t64 \
    libxshmfence1 \
    --no-install-recommends \
    && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------
# Google Chrome
# ------------------------------------------------------------
RUN wget -q https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb \
    && apt-get update \
    && apt-get install -y ./google-chrome-stable_current_amd64.deb \
    && rm google-chrome-stable_current_amd64.deb \
    && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------
# Shiny app
# ------------------------------------------------------------
RUN rm -rf /srv/shiny-server/*
COPY . /srv/shiny-server/

# ------------------------------------------------------------
# R packages
# ------------------------------------------------------------
RUN R -e "install.packages(c( \
  'ggplot2','chromote','dplyr','lubridate','DT','tidyr','plotly', \
  'shiny','kableExtra','reticulate','webshot2','gdata','bslib', \
  'pagedown','shinyMatrix','shinybusy' \
), repos='https://cloud.r-project.org')"

# ------------------------------------------------------------
# Create virtualenv from that Python
# ------------------------------------------------------------
RUN chown -R shiny:shiny /srv/shiny-server
USER shiny
RUN R -e "library(reticulate); \
reticulate::install_python(version = '3.11.9'); \
virtualenv_create('/srv/shiny-server/python_venv'); \
virtualenv_install('/srv/shiny-server/python_venv', \
  packages = c('pydicom==2.4.0', 'rdsr_navigator==0.5.0', 'pandas==2.2.3'))"

# ------------------------------------------------------------
# Port
# ------------------------------------------------------------
EXPOSE 3838

