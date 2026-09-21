$ErrorActionPreference='Stop'
$Root=(Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path; Set-Location $Root
function Invoke-NativeChecked {
    param([string]$FilePath,[string[]]$Arguments,[int[]]$AllowedExitCodes=@(0))
    & $FilePath @Arguments | Out-Host
    $rc=$LASTEXITCODE
    if($AllowedExitCodes -notcontains $rc){throw "native command failed ($rc): $FilePath $($Arguments -join ' ')"}
    return $rc
}
$PythonBin=if($env:PYTHON_BIN){$env:PYTHON_BIN}else{'python'}
Invoke-NativeChecked $PythonBin @('-m','venv','.venv')
$Python=(Join-Path $Root '.venv\Scripts\python.exe')
Invoke-NativeChecked $Python @('-m','pip','install','-e','.[test]')
Invoke-NativeChecked $Python @('-m','proteomics_pipeline','--help')
Invoke-NativeChecked $Python @('-m','proteomics_pipeline','--version')
Invoke-NativeChecked $Python @('-m','proteomics_pipeline','validate','--config','configs/examples/example-independent.json','--schema-only','--json')
$r=Get-Command Rscript -ErrorAction SilentlyContinue
$env:R_LIBS_USER=Join-Path $Root '.r-lib'
$RAvailable=$null -ne $r
if($RAvailable){
    New-Item -ItemType Directory -Force $env:R_LIBS_USER|Out-Null
    Invoke-NativeChecked 'R' @('CMD','INSTALL',"--library=$env:R_LIBS_USER",'r/proteomicsCore')
}
$DoctorAllowed=if($RAvailable){@(0)}else{@(3)}
$DoctorRc=Invoke-NativeChecked $Python @('-m','proteomics_pipeline','doctor','--json') -AllowedExitCodes $DoctorAllowed
Invoke-NativeChecked $Python @('-m','unittest','discover','-s','tests','-v')
Invoke-NativeChecked $Python @('-m','pytest','tests/contract/foundation','tests/unit/foundation','-q')
Invoke-NativeChecked $Python @('scripts/maintained/check_baseline.py')
if(-not $RAvailable){
    Invoke-NativeChecked $Python @('scripts/maintained/write_solved_versions.py') -AllowedExitCodes @(3)
    Write-Output 'NOT_RUN: Rscript is unavailable'
    exit 3
}
Invoke-NativeChecked 'Rscript' @('--vanilla','scripts/maintained/test_r.R','--suite','foundation')
Invoke-NativeChecked $Python @('scripts/maintained/write_solved_versions.py')
exit $DoctorRc
