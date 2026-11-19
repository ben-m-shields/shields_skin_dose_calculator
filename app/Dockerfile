FROM rocker/shiny-verse:latest

# ------------------------------------------------------------
# System dependencies
# ------------------------------------------------------------
# GOOGLE CHROME AND ITS FOURTEEN X/GTK LIBRARIES ARE GONE, and nothing replaced them. They existed for one
# reason: report.Rmd rendered Figure 2 by saving a plotly mesh3d as an HTML widget, booting Chrome, letting
# plotly.js build a WebGL context (software-rasterised in a container, which is the slow part) and
# screenshotting the page. skin_dose_map.R now rasterises the same mesh in R with a painter's algorithm -
# measured at ~1 s for 82,000 triangles across four views, against a browser launch per report - so the
# libcups2 / libgtk-3-0 / libxcomposite1 / libxfixes3 / xdg-utils / libnss3 / libatk-bridge2.0-0 /
# libxkbcommon0 / libxdamage1 / libxrandr2 / libgbm1 / libasound2t64 / libxshmfence1 / fonts-liberation
# set and the ~400 MB Chrome .deb all came out with it.
#
# Do NOT reintroduce a browser (or a LaTeX/pandoc toolchain) to make a document. officer writes .docx with
# no external toolchain, and base R rasterises a triangle mesh perfectly well.
RUN apt-get update && apt-get install -y \
	patch \
    curl \
    git \
    wget \
    --no-install-recommends \
    && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------
# Shiny app
# ------------------------------------------------------------
RUN rm -rf /srv/shiny-server/*
COPY . /srv/shiny-server/

# ------------------------------------------------------------
# R packages
# ------------------------------------------------------------
# chromote / webshot2 / pagedown / kableExtra left with report.Rmd (see the Chrome note above);
# officer + flextable build the Word report, ragg rasterises Figure 2. rocker/shiny-verse already carries
# rmarkdown and knitr, but nothing in the app calls them any more.
RUN R -e "install.packages(c( \
  'ggplot2','dplyr','lubridate','DT','tidyr','plotly', \
  'shiny','reticulate','gdata','bslib', \
  'shinyMatrix','shinybusy','officer','flextable','ragg' \
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

