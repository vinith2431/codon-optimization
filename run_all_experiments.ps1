$ErrorActionPreference = "Continue"

$logDir = "results\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$logFile = Join-Path $logDir ("full_run_" + (Get-Date -Format "yyyyMMdd_HHmmss") + ".log")

function Run-Step {
    param([string]$Title, [string]$Command)

    Write-Host "`n==================================================" -ForegroundColor Cyan
    Write-Host $Title -ForegroundColor Cyan
    Write-Host "==================================================" -ForegroundColor Cyan

    "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] START: $Title" |
        Tee-Object -FilePath $logFile -Append

    Invoke-Expression "$Command 2>&1" |
        Tee-Object -FilePath $logFile -Append

    $code = $LASTEXITCODE
    if ($null -eq $code) { $code = 0 }

    if ($code -eq 0) {
        Write-Host "PASSED: $Title" -ForegroundColor Green
        "RESULT: PASS" | Tee-Object -FilePath $logFile -Append
    } else {
        Write-Host "FAILED (exit $code): $Title" -ForegroundColor Red
        "RESULT: FAIL (exit $code)" | Tee-Object -FilePath $logFile -Append
    }
}

Write-Host "Python environment:" -ForegroundColor Yellow
python --version
python -c "import sys; print(sys.executable)"

Run-Step "1. Unit tests" "python -m pytest -q"
Run-Step "2. Full optimizer comparison" "python experiments\run_optimizer_comparison.py"
Run-Step "3. Research experiment suite" "python experiments\run_research_suite.py"
Run-Step "4. Constraint-aware multi-start benchmark" "python experiments\run_constraint_aware_multistart.py"
Run-Step "5. Pareto benchmark" "python experiments\run_pareto_benchmark.py"

Write-Host "`n==================================================" -ForegroundColor Cyan
Write-Host "GPU / CUDA CHECK" -ForegroundColor Cyan
Write-Host "==================================================" -ForegroundColor Cyan

python -c "import torch; print('PyTorch:', torch.__version__); print('CUDA available:', torch.cuda.is_available()); print('GPU:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'Not available')" 2>&1 |
    Tee-Object -FilePath $logFile -Append

Write-Host "`nAll steps attempted." -ForegroundColor Green
Write-Host "Log file: $logFile" -ForegroundColor Yellow
Write-Host "Review any FAILED steps before trusting the combined results." -ForegroundColor Yellow
