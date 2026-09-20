# Configure Python Virtual Environment for FDS (WINDOWS)
# Usage (PowerShell):
#    . .\setup_python_env.ps1 [--batchmode]
#
# Dot-source this script (leading ". ") so activation is guaranteed to remain
# in the current PowerShell session.

function Invoke-FdsPythonEnvSetup {
    param(
        [bool]$BatchMode
    )

    $OriginalDirectory = (Get-Location).Path

    try {
        # === Check for python ===
        if (-not (Get-Command python -ErrorAction SilentlyContinue)) {
            throw "python is not installed or not in PATH"
        }

        # === Ensure Python version >= 3.7 ===
        $VersionText = & python -c "import sys; print(sys.version_info.major, sys.version_info.minor)"
        if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($VersionText)) {
            throw "Failed to determine Python version"
        }

        $VersionParts = $VersionText.Trim().Split([char[]]' ', [System.StringSplitOptions]::RemoveEmptyEntries)
        if ($VersionParts.Count -lt 2) {
            throw "Failed to parse Python version: $VersionText"
        }

        $PyMajor = [int]$VersionParts[0]
        $PyMinor = [int]$VersionParts[1]

        if ($PyMajor -lt 3 -or ($PyMajor -eq 3 -and $PyMinor -lt 7)) {
            throw "Python 3.7 or higher is required. Found $PyMajor.$PyMinor."
        }

        Write-Host "Python version is OK: $PyMajor.$PyMinor"

        # === Locate repo root and .github folder ===
        $RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..\..")).Path
        $GitHubDir = Join-Path $RepoRoot "fds\.github"
        if (-not (Test-Path -LiteralPath $GitHubDir -PathType Container)) {
            throw "Directory not found: $GitHubDir"
        }

        Set-Location -LiteralPath $GitHubDir

        # === Setup virtual environment ===
        $VenvDir = Join-Path $GitHubDir "fds_python_env"
        $InstallRequirements = $false
        $Reinstall = $false

        if (Test-Path -LiteralPath $VenvDir -PathType Container) {
            if ($BatchMode) {
                # Preserve existing Windows batch-mode behavior: do not delete an existing environment.
                Write-Host "Batch mode: activating existing virtual environment without prompts or deletion."
            }
            else {
                Write-Host "Virtual environment '$VenvDir' already exists."
                $Choice = Read-Host "Do you want to reinstall everything? (y/N)"
                if ($Choice -match '^(?i:y|yes)$') {
                    $Reinstall = $true
                }
                else {
                    Write-Host "Activating existing environment..."
                }
            }
        }
        else {
            $Reinstall = $true
        }

        if ($Reinstall) {
            # Deactivate any currently active venv before removing/recreating this one.
            if ($env:VIRTUAL_ENV) {
                $DeactivateCommand = Get-Command deactivate -ErrorAction SilentlyContinue
                if ($DeactivateCommand) {
                    deactivate
                }
            }

            if (Test-Path -LiteralPath $VenvDir -PathType Container) {
                Write-Host "Removing old environment..."
                Remove-Item -LiteralPath $VenvDir -Recurse -Force -ErrorAction Stop
            }

            Write-Host "Creating new virtual environment..."
            & python -m venv $VenvDir
            if ($LASTEXITCODE -ne 0) {
                throw "Failed to create virtual environment"
            }
            $InstallRequirements = $true
        }

        # === Activate environment ===
        $ActivateScript = Join-Path $VenvDir "Scripts\Activate.ps1"
        if (-not (Test-Path -LiteralPath $ActivateScript -PathType Leaf)) {
            throw "Virtual environment activation script not found: $ActivateScript"
        }

        # Dot-source Activate.ps1 so the environment remains active in this shell.
        . $ActivateScript

        if (-not $env:VIRTUAL_ENV) {
            throw "Failed to activate virtual environment"
        }

        # === Upgrade pip and install requirements if flagged ===
        if ($InstallRequirements) {
            Write-Host "Installing/updating required Python packages..."
            & python -m pip install --upgrade pip
            if ($LASTEXITCODE -ne 0) {
                throw "Failed to upgrade pip"
            }

            $RequirementsFile = Join-Path $GitHubDir "requirements.txt"
            if (Test-Path -LiteralPath $RequirementsFile -PathType Leaf) {
                & python -m pip install -r $RequirementsFile
                if ($LASTEXITCODE -ne 0) {
                    throw "Failed to install requirements"
                }
            }
        }

        # === Set PYTHONPATH ===
        $FdsPythonPath = Join-Path $RepoRoot "fds\Utilities\Python"
        if ([string]::IsNullOrWhiteSpace($env:PYTHONPATH)) {
            $env:PYTHONPATH = $FdsPythonPath
        }
        else {
            $env:PYTHONPATH = "$FdsPythonPath;$env:PYTHONPATH"
        }

        # === Run test script ===
        if (-not (Test-Path -LiteralPath $FdsPythonPath -PathType Container)) {
            throw "Failed to find script directory: $FdsPythonPath"
        }

        Set-Location -LiteralPath $FdsPythonPath
        if (-not (Test-Path -LiteralPath "hello_world.py" -PathType Leaf)) {
            throw "hello_world.py not found"
        }

        & python .\hello_world.py
        if ($LASTEXITCODE -ne 0) {
            throw "hello_world.py failed"
        }

        Write-Host ""
        Write-Host "Python environment setup complete."
    }
    finally {
        Set-Location -LiteralPath $OriginalDirectory
    }
}

# Accept both the existing --batchmode spelling and the PowerShell-style -BatchMode.
$BatchModeRequested = $false
foreach ($Argument in $args) {
    if ($Argument -ieq "--batchmode" -or $Argument -ieq "-BatchMode") {
        $BatchModeRequested = $true
    }
    else {
        throw "Unknown argument: $Argument"
    }
}

try {
    Invoke-FdsPythonEnvSetup -BatchMode:$BatchModeRequested
}
catch {
    throw "*** Error: $($_.Exception.Message)"
}
finally {
    Remove-Item Function:\Invoke-FdsPythonEnvSetup -ErrorAction SilentlyContinue
    Remove-Variable BatchModeRequested -ErrorAction SilentlyContinue
    Remove-Variable Argument -ErrorAction SilentlyContinue
}
