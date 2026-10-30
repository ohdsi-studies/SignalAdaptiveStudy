# Package repository: Posit Package Manager binaries. This file replaces ~/.Rprofile inside the project, so the
# repository is set here; without it, installs and renv::restore() fall back to compiling CRAN source packages.
# Linux: P3M's generic manylinux build (as ~/.Rprofile); elsewhere P3M's standard URL serves Windows/macOS binaries.
# Set before activation for renv (repos.override) and again after it, because renv resets repos from renv.lock.
local({
  p3m <- if (Sys.info()[["sysname"]] == "Linux") "https://p3m.dev/cran/__linux__/manylinux_2_28/latest"
         else "https://p3m.dev/cran/latest"
  options(renv.config.repos.override = p3m)
  source("renv/activate.R")
  options(repos = c(CRAN = p3m))
})
