# Pre-Push Verification and Git Commands
# Run this before pushing to GitHub

Write-Host "═══════════════════════════════════════════════════════" -ForegroundColor Cyan
Write-Host "  BTEH: Pre-Push Verification" -ForegroundColor Cyan
Write-Host "═══════════════════════════════════════════════════════" -ForegroundColor Cyan
Write-Host ""

# Navigate to project root
Set-Location "d:\PhD\Projects\SDM_projects\BTEH\SDM_ele\BTEH"

Write-Host "✓ Current directory: $(Get-Location)" -ForegroundColor Green
Write-Host ""

# Stage all new documentation
Write-Host "Adding documentation files..." -ForegroundColor Yellow
git add *.md
git add *.R
git add *.sh

# Stage lock files
Write-Host "Adding renv lock files..." -ForegroundColor Yellow
git add renv.lock renv_*.lock

# Stage modified scripts
Write-Host "Adding modified scripts..." -ForegroundColor Yellow
git add "BTEH/scripts_RemoteSensing for E and C_/*.R"
git add BTEH/scripts_ems_paper/*.R

# Stage other changes
git add README.md

Write-Host ""
Write-Host "═══════════════════════════════════════════════════════" -ForegroundColor Cyan
Write-Host "  Files Staged for Commit" -ForegroundColor Cyan
Write-Host "═══════════════════════════════════════════════════════" -ForegroundColor Cyan
git status --short

Write-Host ""
Write-Host "═══════════════════════════════════════════════════════" -ForegroundColor Cyan
Write-Host "  Verification Checks" -ForegroundColor Cyan
Write-Host "═══════════════════════════════════════════════════════" -ForegroundColor Cyan

# Check 1: Lock file
Write-Host ""
Write-Host "Check 1: Verify renv.lock is Remote Sensing version" -ForegroundColor Yellow
$lockSize = (Get-Item renv.lock).Length
$rsLockSize = (Get-Item renv_remotesensing.lock).Length
if ($lockSize -eq $rsLockSize) {
    Write-Host "  ✓ renv.lock matches renv_remotesensing.lock" -ForegroundColor Green
} else {
    Write-Host "  ⚠ renv.lock size: $lockSize" -ForegroundColor Yellow
    Write-Host "  ⚠ renv_remotesensing.lock size: $rsLockSize" -ForegroundColor Yellow
}

# Check 2: Documentation
Write-Host ""
Write-Host "Check 2: Documentation files present" -ForegroundColor Yellow
$docs = @(
    "DUAL_STUDY_GUIDE.md",
    "END_TO_END_VERIFICATION.md",
    "HPC_DEPLOYMENT_CHECKLIST.md",
    "GITHUB_HPC_DEPLOYMENT.md"
)
foreach ($doc in $docs) {
    if (Test-Path $doc) {
        Write-Host "  ✓ $doc" -ForegroundColor Green
    } else {
        Write-Host "  ✗ $doc MISSING" -ForegroundColor Red
    }
}

# Check 3: Helper scripts
Write-Host ""
Write-Host "Check 3: Helper scripts present" -ForegroundColor Yellow
$scripts = @(
    "create_study_locks.R",
    "setup_for_hpc.R",
    "run_pipeline_hpc.sh"
)
foreach ($script in $scripts) {
    if (Test-Path $script) {
        Write-Host "  ✓ $script" -ForegroundColor Green
    } else {
        Write-Host "  ✗ $script MISSING" -ForegroundColor Red
    }
}

# Check 4: No large files
Write-Host ""
Write-Host "Check 4: Checking for large files (>10MB)..." -ForegroundColor Yellow
$largeFiles = git ls-files | ForEach-Object {
    if (Test-Path $_) {
        $size = (Get-Item $_).Length
        if ($size -gt 10MB) {
            [PSCustomObject]@{
                File = $_
                Size = "{0:N2} MB" -f ($size / 1MB)
            }
        }
    }
}
if ($largeFiles) {
    Write-Host "  ⚠ Large files found:" -ForegroundColor Yellow
    $largeFiles | Format-Table -AutoSize
} else {
    Write-Host "  ✓ No large files" -ForegroundColor Green
}

# Check 5: renv/library not included
Write-Host ""
Write-Host "Check 5: Verify renv/library is not staged" -ForegroundColor Yellow
$staged = git diff --cached --name-only | Select-String "renv/library"
if ($staged) {
    Write-Host "  ✗ ERROR: renv/library is staged!" -ForegroundColor Red
    Write-Host "  Run: git reset renv/library" -ForegroundColor Yellow
} else {
    Write-Host "  ✓ renv/library not staged" -ForegroundColor Green
}

Write-Host ""
Write-Host "═══════════════════════════════════════════════════════" -ForegroundColor Cyan
Write-Host "  Ready to Commit?" -ForegroundColor Cyan
Write-Host "═══════════════════════════════════════════════════════" -ForegroundColor Cyan
Write-Host ""
Write-Host "Next steps:" -ForegroundColor Yellow
Write-Host "  1. Review staged files above" -ForegroundColor White
Write-Host "  2. Run commit command (see below)" -ForegroundColor White
Write-Host "  3. Push to GitHub" -ForegroundColor White
Write-Host "  4. Pull on HPC" -ForegroundColor White
Write-Host ""
Write-Host "Commit command:" -ForegroundColor Green
Write-Host '  git commit -m "HPC deployment ready: Remote Sensing (123 packages)"' -ForegroundColor Cyan
Write-Host ""
Write-Host "Push command:" -ForegroundColor Green
Write-Host "  git push origin main" -ForegroundColor Cyan
Write-Host ""
