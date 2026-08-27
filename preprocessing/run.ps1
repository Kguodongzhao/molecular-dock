[CmdletBinding()]
param(
    [Parameter(Mandatory = $true, Position = 0)]
    [ValidateSet("ligands", "receptor")]
    [string]$Action,

    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$RemainingArguments
)

$ErrorActionPreference = "Stop"
$scriptRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$pythonScript = if ($Action -eq "ligands") {
    Join-Path $scriptRoot "prepare_ligands.py"
} else {
    Join-Path $scriptRoot "prepare_receptor.py"
}

$condaCommand = Get-Command conda -ErrorAction SilentlyContinue
$condaPath = if ($condaCommand) { $condaCommand.Source } else { $null }
if (-not $condaPath) {
    $candidates = @(
        "C:\ProgramData\anaconda3\condabin\conda.bat",
        (Join-Path $env:USERPROFILE "anaconda3\condabin\conda.bat"),
        (Join-Path $env:USERPROFILE "miniconda3\condabin\conda.bat"),
        (Join-Path $env:USERPROFILE "miniforge3\condabin\conda.bat")
    )
    $condaPath = $candidates | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
}
if (-not $condaPath) {
    throw "Cannot find conda. Activate the molecular-docking environment and run the Python script directly."
}

Write-Host "Using conda environment: molecular-docking"
& $condaPath run --no-capture-output -n molecular-docking python $pythonScript @RemainingArguments
exit $LASTEXITCODE
