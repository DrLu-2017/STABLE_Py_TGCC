# run-claude-ollama.ps1
# 远程 Ollama 服务版本，不需要本地安装 ollama

#Requires -Version 5.1
$ErrorActionPreference = "Stop"

function Write-Ok   { param($msg) Write-Host "  [OK] $msg" -ForegroundColor Green }
function Write-Warn { param($msg) Write-Host "  [!!] $msg" -ForegroundColor Yellow }
function Write-Title{ param($msg) Write-Host "`n$msg" -ForegroundColor White }

$Config = @{
    OllamaBaseUrl   = "http://127.0.0.1:11434"
    OllamaAuthToken = "ollama"
    DefaultModel    = "qwen3-coder:30b"       # ← 更新为实际可用模型
    HaikuModel      = "gemma4-26b-custom:latest"       # ← 轻量任务（原 Haiku 角色）
    SonnetModel     = "qwen2.5-coder-32k:latest"       # ← 标准任务（原 Sonnet 角色）
    OpusModel       = "glm-4.7-flash:q8_0"  # ← 重型任务（原 Opus 角色，31GB）
}

Clear-Host
Write-Host "==================================================" -ForegroundColor DarkCyan
Write-Host "   Claude Code + Remote Ollama Launch Script     " -ForegroundColor DarkCyan
Write-Host "==================================================" -ForegroundColor DarkCyan

# 1. Check remote Ollama service
Write-Title "[ 1/3 ] Checking remote Ollama at $($Config.OllamaBaseUrl)"

try {
    $resp = Invoke-RestMethod -Uri "$($Config.OllamaBaseUrl)/api/tags" -TimeoutSec 5 -ErrorAction Stop
    Write-Ok "Service reachable, $($resp.models.Count) model(s) available"

    $installedModels = $resp.models | ForEach-Object { $_.name }

    # 显示所有可用模型
    Write-Host "`n  Available models:" -ForegroundColor Cyan
    $installedModels | ForEach-Object { Write-Host "       - $_" -ForegroundColor DarkGray }

    if ($installedModels -contains $Config.DefaultModel) {
        Write-Ok "Default model ready: $($Config.DefaultModel)"
    } else {
        Write-Warn "Model '$($Config.DefaultModel)' not found. Available models listed above."
        Read-Host "Press Enter to exit"
        exit
    }
} catch {
    Write-Warn "Cannot reach Ollama at $($Config.OllamaBaseUrl)"
    Read-Host "Press Enter to exit"
    exit
}

# 2. Check Claude CLI
Write-Title "[ 2/3 ] Checking Claude CLI"

$claude = Get-Command claude -ErrorAction SilentlyContinue
if (-not $claude) {
    Write-Warn "claude not found. Install with: npm install -g @anthropic-ai/claude-code"
    Read-Host "Press Enter to exit"
    exit
}
Write-Ok "Claude CLI: $($claude.Source)"

# 3. Set env vars and launch
Write-Title "[ 3/3 ] Setting environment variables and launching"
# 设置环境变量（限制上限调大至 64000）
$env:CLAUDE_CODE_MAX_OUTPUT_TOKENS = "64000"
# 如果中间件或服务支持透传上下文设置
$env:OLLAMA_NUM_PARALLEL = "1"

$env:ANTHROPIC_BASE_URL             = $Config.OllamaBaseUrl
$env:ANTHROPIC_AUTH_TOKEN           = $Config.OllamaAuthToken
$env:ANTHROPIC_API_KEY              = "dummy"
$env:ANTHROPIC_MODEL                = $Config.DefaultModel
$env:ANTHROPIC_DEFAULT_HAIKU_MODEL  = $Config.HaikuModel
$env:ANTHROPIC_DEFAULT_SONNET_MODEL = $Config.SonnetModel
$env:ANTHROPIC_DEFAULT_OPUS_MODEL   = $Config.OpusModel
$env:CLAUDE_CODE_USE_BEDROCK        = "0"
$env:CLAUDE_CODE_USE_VERTEX         = "0"

Write-Ok "Default : $($Config.DefaultModel)"
Write-Ok "Haiku   : $($Config.HaikuModel)"
Write-Ok "Sonnet  : $($Config.SonnetModel)"
Write-Ok "Opus    : $($Config.OpusModel)"
Write-Ok "URL     : $($Config.OllamaBaseUrl)"
Write-Host ""

claude

Write-Host ""
Write-Ok "Session ended."
Read-Host "Press Enter to close"