$ErrorActionPreference = "Stop"

$RootDir = Resolve-Path (Join-Path $PSScriptRoot "..")
Set-Location $RootDir

conda env create -f environment.yml

Write-Host ""
Write-Host "Conda environment created."
Write-Host ""
Write-Host "Next:"
Write-Host "  conda activate pymol-siteview"
Write-Host "  pymol-siteview --help"
Write-Host ""
Write-Host "Example:"
Write-Host "  pymol-siteview your_structure.pdb -o siteview_out"
