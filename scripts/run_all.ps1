# Run the whole study from step 1 to step 7, then the tests.
# Usage (from the repository root):  powershell -File scripts\run_all.ps1
$ErrorActionPreference = "Stop"
$py = ".\.venv\Scripts\python.exe"
foreach ($s in "01_build_panels", "02_screen_guides", "03_independent_check", "04_make_figure",
               "05_genome_pilot", "06_allele_misses", "07_report_tables") {
    Write-Host "=== $s ==="
    & $py "scripts\$s.py"
    if ($LASTEXITCODE -ne 0) { throw "$s failed with exit code $LASTEXITCODE" }
}
Write-Host "=== pytest ==="
& $py -m pytest -q
if ($LASTEXITCODE -ne 0) { throw "tests failed" }
