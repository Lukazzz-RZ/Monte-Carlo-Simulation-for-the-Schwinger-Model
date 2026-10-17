# Bounded execution check. The original source files are NEVER modified:
# a temporary copy is made and the run constants are shortened there.
# This checks execution, not the physical validity of the algorithm.
$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$root = Join-Path ([IO.Path]::GetTempPath()) ("schwinger-smoke-" + [guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Path $root -Force | Out-Null

try {
    Copy-Item (Join-Path $repoRoot "CMakeLists.txt") $root
    Copy-Item (Join-Path $repoRoot "include") $root -Recurse
    Copy-Item (Join-Path $repoRoot "src") $root -Recurse

    $mainPath = Join-Path $root "src/Main.c"
    $source = [IO.File]::ReadAllText($mainPath)

    $replacements = @{
        "static const long N_THERMALIZATION = 100000L;" = "static const long N_THERMALIZATION = 4L;"
        "static const long N_MEASUREMENTS = 50176L;" = "static const long N_MEASUREMENTS = 256L;"
        "static const int STEPS_BETWEEN_MEASUREMENTS = 10;" = "static const int STEPS_BETWEEN_MEASUREMENTS = 1;"
        "static const long PROGRESS_EVERY = 4096L;" = "static const long PROGRESS_EVERY = 256L;"
    }
    foreach ($old in $replacements.Keys) {
        if (-not $source.Contains($old)) {
            throw "Expected setting not found in src/Main.c: $old"
        }
        $source = $source.Replace($old, $replacements[$old])
    }

    $thetaPattern = 'static const double THETA_OVER_PI_VALUES\[\]\s*=\s*\{[\s\S]*?\};'
    if ([regex]::Matches($source, $thetaPattern).Count -ne 1) {
        throw "Could not safely identify the theta sweep in src/Main.c"
    }
    $source = [regex]::Replace(
        $source, $thetaPattern,
        'static const double THETA_OVER_PI_VALUES[] = { 0.0 };'
    )
    [IO.File]::WriteAllText($mainPath, $source, (New-Object Text.UTF8Encoding($false)))

    & cmake -S $root -B (Join-Path $root "build") -G "Visual Studio 17 2022" -A x64
    if ($LASTEXITCODE -ne 0) { throw "Smoke CMake configuration failed" }

    & cmake --build (Join-Path $root "build") --config Release
    if ($LASTEXITCODE -ne 0) { throw "Smoke CMake build failed" }

    $exe = Join-Path $root "build/Release/SchwingerWorm.exe"
    if (-not (Test-Path $exe)) { throw "Smoke executable missing" }
    $out = Join-Path $root "run.stdout.txt"
    $err = Join-Path $root "run.stderr.txt"

    $process = Start-Process -FilePath $exe -WorkingDirectory $root -RedirectStandardOutput $out -RedirectStandardError $err -PassThru
    if (-not $process.WaitForExit(90000)) {
        $process.Kill($true)
        throw "Smoke simulation exceeded 90 seconds"
    }
    $process.Refresh()
    Get-Content $out | Select-Object -Last 12 | Write-Host
    if ((Test-Path $err) -and (Get-Item $err).Length -gt 0) {
        Get-Content $err | Select-Object -Last 30 | Write-Host
    }
    if ($process.ExitCode -ne 0) {
        throw "Simulation exited with code $($process.ExitCode)"
    }

    $resultDir = Join-Path $root "Results"
    $series = @(Get-ChildItem -Path $resultDir -Filter "Series_*.dat" -File)
    if ($series.Count -ne 1) {
        throw "Expected exactly one measurement file; found $($series.Count)"
    }
    $rows = @(Get-Content $series[0].FullName | Where-Object { $_ -and -not $_.StartsWith("#") })
    if ($rows.Count -ne 256) { throw "Expected 256 measurements; found $($rows.Count)" }
    foreach ($row in $rows) {
        $values = @($row.Trim() -split '\s+')
        if ($values.Count -ne 22) { throw "Incorrect column count: $($values.Count)" }
        if ($row -match '(?i)(nan|inf)') { throw "Non-finite data in measurement file" }
    }

    $thermal = @(Get-ChildItem -Path $resultDir -Filter "Thermalization_*.txt" -File)
    if ($thermal.Count -ne 1 -or $thermal[0].Length -eq 0) {
        throw "Missing or empty thermalisation record"
    }
    Write-Host "Smoke test passed: 256 finite measurements, 22 columns each."
}
finally {
    Remove-Item -LiteralPath $root -Recurse -Force -ErrorAction SilentlyContinue
}

