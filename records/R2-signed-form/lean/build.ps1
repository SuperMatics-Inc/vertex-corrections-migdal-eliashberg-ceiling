$ErrorActionPreference = 'Stop'
Push-Location $PSScriptRoot
try {
    $packages = 'C:/dev/.cache/lean/mathlib-v4.35.0/.lake/packages'
    $overlay = Join-Path $PSScriptRoot 'prebuilt/.lake/build/lib/lean'
    New-Item -ItemType Directory -Force $overlay | Out-Null
    foreach ($package in (Get-ChildItem -Directory $packages)) {
        $library = Join-Path $package.FullName '.lake/build/lib/lean'
        if (!(Test-Path -LiteralPath $library)) { continue }
        foreach ($entry in (Get-ChildItem -LiteralPath $library)) {
            $destination = Join-Path $overlay $entry.Name
            if ($entry.PSIsContainer) {
                if (!(Test-Path -LiteralPath $destination)) {
                    New-Item -ItemType Junction -Path $destination -Target $entry.FullName | Out-Null
                }
            } elseif ($entry.Name -match '\.(olean|ir)(\.(private|server|sig))?$') {
                Copy-Item -LiteralPath $entry.FullName -Destination $destination
            }
        }
    }
    if (!(Test-Path -LiteralPath (Join-Path $overlay 'Mathlib/Analysis/Matrix/PosDef.olean'))) {
        throw 'Missing the shared Mathlib v4.35.0 artifacts.'
    }
    & lake build
    exit $LASTEXITCODE
} finally {
    Pop-Location
}
