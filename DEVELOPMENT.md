# Development notes

## First local setup

Open `suest.Rproj`, then run:

```r
install.packages(c(
  "devtools",
  "pkgdown",
  "roxygen2",
  "testthat"
))

devtools::document()
devtools::test()
devtools::check()
pkgdown::build_site()
```

## Create the GitHub repository

Create a new public repository named `suest` under `tdmize`, then from the
package directory run:

```bash
git init
git add .
git commit -m "Initial suest package"
git branch -M main
git remote add origin https://github.com/tdmize/suest.git
git push -u origin main
```

The included GitHub Actions workflows will:

- run `R CMD check` on Windows, macOS, and Linux;
- build and deploy the pkgdown site to the `gh-pages` branch;
- run the full Mize, Doan, and Long (2019) numerical acceptance suite.

After the first successful pkgdown deployment, enable GitHub Pages using the
`gh-pages` branch if GitHub has not done so automatically.
