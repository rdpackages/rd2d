#!/usr/bin/env Rscript

options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
mode <- "release"
if ("--dev" %in% args) mode <- "dev"
if ("--pre-push" %in% args) mode <- "pre-push"
if ("--release" %in% args) mode <- "release"
keep_artifacts <- "--keep-artifacts" %in% args

find_repo_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "R", "rd2d", "DESCRIPTION"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not find rd2d repository root.", call. = FALSE)
    path <- parent
  }
}

repo_root <- find_repo_root()
pkg_dir <- file.path(repo_root, "R", "rd2d")

add_common_paths <- function() {
  if (.Platform$OS.type == "windows") {
    Sys.unsetenv(c("LC_ALL", "LC_CTYPE"))
    r_minor <- sub("\\..*", "", R.version$minor)
    local_app_data <- chartr("\\", "/", Sys.getenv("LOCALAPPDATA"))
    user_lib <- file.path(local_app_data, "R", "win-library",
                          paste(R.version$major, r_minor, sep = "."))
    if (dir.exists(user_lib)) {
      Sys.setenv(R_LIBS_USER = user_lib)
      Sys.setenv(R_LIBS = user_lib)
      .libPaths(c(user_lib, .libPaths()))
    }
  }
  if (.Platform$OS.type == "windows" && !nzchar(Sys.which("pandoc"))) {
    pandoc_dir <- file.path(Sys.getenv("LOCALAPPDATA"), "Pandoc")
    pandoc_exe <- file.path(pandoc_dir, "pandoc.exe")
    if (file.exists(pandoc_exe)) {
      Sys.setenv(PATH = paste(pandoc_dir, Sys.getenv("PATH"), sep = .Platform$path.sep))
    }
  }
}

run <- function(command, args = character(), wd = repo_root) {
  cat("\n==>", command, paste(args, collapse = " "), "\n")
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(wd)
  child_env <- character()
  if (nzchar(Sys.getenv("R_LIBS_USER"))) {
    child_env <- c(child_env, paste0("R_LIBS_USER=", Sys.getenv("R_LIBS_USER")))
  }
  if (nzchar(Sys.getenv("R_LIBS"))) {
    child_env <- c(child_env, paste0("R_LIBS=", Sys.getenv("R_LIBS")))
  }
  status <- system2(command, args, env = child_env)
  if (!identical(status, 0L)) {
    stop(sprintf("Command failed with status %s: %s", status, command), call. = FALSE)
  }
  invisible(TRUE)
}

run_r <- function(expr, wd = repo_root) {
  run(file.path(R.home("bin"), "Rscript"), c("-e", expr), wd = wd)
}

require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing required package(s): ", paste(missing, collapse = ", "),
      "\nInstall them in your user library before running this check.",
      call. = FALSE
    )
  }
}

clean_artifacts <- function() {
  for (attempt in seq_len(5)) {
    paths <- c(
      file.path(repo_root, "check_as_cran"),
      file.path(repo_root, "rd2d.Rcheck"),
      Sys.glob(file.path(repo_root, "rd2d_*.tar.gz"))
    )
    paths <- paths[file.exists(paths)]
    if (!length(paths)) return(invisible(TRUE))
    unlink(paths, recursive = TRUE, force = TRUE)
    if (!any(file.exists(paths))) return(invisible(TRUE))
    Sys.sleep(1)
  }
  remaining <- paths[file.exists(paths)]
  if (length(remaining)) {
    message("Could not remove generated artifact(s): ", paste(remaining, collapse = ", "))
  }
  invisible(FALSE)
}

add_common_paths()

cat("Repository:", repo_root, "\n")
cat("Package:   ", pkg_dir, "\n")
cat("Mode:      ", mode, "\n")
cat("Pandoc:    ", if (nzchar(Sys.which("pandoc"))) Sys.which("pandoc") else "<not found>", "\n")

if (mode == "pre-push") {
  run(file.path(R.home("bin"), "R"), c("CMD", "check", "--no-manual", "R/rd2d"))
  if (!keep_artifacts) clean_artifacts()
  cat("\nPre-push checks passed.\n")
  quit(status = 0)
}

require_packages(c("roxygen2", "testthat"))
run_r("roxygen2::roxygenise()", wd = pkg_dir)
run_r("testthat::test_local()", wd = pkg_dir)

if (mode == "dev") {
  run(file.path(R.home("bin"), "R"), c("CMD", "check", "--no-manual", "R/rd2d"))
  if (!keep_artifacts) clean_artifacts()
  cat("\nDevelopment checks passed.\n")
  quit(status = 0)
}

run(file.path(R.home("bin"), "R"), c("CMD", "build", "R/rd2d"))
tarball <- Sys.glob(file.path(repo_root, "rd2d_*.tar.gz"))
if (length(tarball) != 1L) stop("Expected exactly one rd2d source tarball.", call. = FALSE)

dir.create(file.path(repo_root, "check_as_cran"), showWarnings = FALSE)
run(
  file.path(R.home("bin"), "R"),
  c("CMD", "check", "--no-manual", "--as-cran", "-o", "check_as_cran", basename(tarball))
)

if (!keep_artifacts) clean_artifacts()
cat("\nRelease checks passed.\n")
