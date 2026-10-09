# 公共环境的版本管理

这里是 16 个 scRNA skill 使用的六组环境配置的规范来源。公共运行目录默认是 `~/projects/scrna_envs`；其中 `.pixi/` 是本机安装结果，不作为源码提交。仅支持锁文件中的 `linux-64` 平台。

| 环境项目 | Skill | 部署的子环境 |
|---|---|---|
| 01-scrna-qc | 01–04 | default |
| 02-annotation | 07–09、14–15 | default |
| 03-integration | 05–06 | default、scvi |
| 05-pathway_program | 10、16 | default、cnmf |
| 06-deg-analysis | 11–12 | default |
| 07-cell-abundance | 13 | default、sccoda |

`pixi.toml` 声明依赖，`pixi.lock` 固定解析结果；整合项目还保存镜像 channel 的 Conda/PyPI 映射。`bundle.json` 记录配置 SHA256 和需要部署的子环境。评分项目中的 cnmf 子环境供 16 的程序发现使用，部署工具会安装它；decoupler 保留原声明和锁定记录，当前执行器不使用，部署工具不会默认安装它。

## 核查与部署

在仓库根目录执行：

```bash
# 只读核查，不创建目录、不安装依赖
python3 scripts/manage_environments.py --target ~/projects/scrna_envs

# 新机器：复制配置并显式安装依赖
python3 scripts/manage_environments.py --target ~/projects/scrna_envs --apply --install

# 已有机器：先备份配置，再同步；不修改二进制环境
python3 scripts/manage_environments.py --target ~/projects/scrna_envs --apply --force

# 单独部署丰度环境
python3 scripts/manage_environments.py --target ~/projects/scrna_envs \
  --profile 07-cell-abundance --apply --force --install
```

备份放在目标根目录的 `.environment-backups/`。工具先验证全部源配置的校验值，再修改目标；复制失败会恢复已替换的文件。依赖安装失败时保留配置备份并返回错误，不承诺回滚整个二进制环境。需要隔离旧环境时，部署到新目录并通过 `runtime.pixi_root` 指向新目录。

`--install` 调用 `pixi install --frozen --run-post-link-scripts`，让锁定的 Bioconductor 数据包完成安装步骤；随后在已激活的 R 环境运行补装程序。这些操作仅用于显式部署。配置审计不能证明所有 Conda/PyPI 包未被人工修改；选定方法仍须运行 skill 的依赖检查和实际分析验证。

## 补装 R 包

`supplemental-r.json` 固定准确版本、源归档和校验值：

- Bioconductor 注释数据使用现有 Bioconda `dataURLs.json` 的版本和 MD5，避免遗漏安装后的下载步骤；评分环境中的 Reactome 数据也纳入清单。
- 丰度环境的 checkmate、distributional、tensorA 使用实测源归档 SHA256。
- posterior、cmdstanr 固定完整 Git 提交和实测归档 SHA256。

`install_supplemental.R --check` 只核查目标环境库中的版本及已记录的 Git 来源；`--install` 才下载、校验和安装。正确版本会跳过。源码安装使用 `dependencies=FALSE`，依赖由锁文件和清单顺序提供，避免滚动升级。Git 归档的安装来源写入包目录的 `.scrna-source.json`。安装时通过 `pixi run --frozen --no-install ... -e default -- Rscript ...` 激活编译器和构建环境。部署工具已经包含这一步。

补装源归档未包含在 Git 仓库中。换机器仍需要下载 Conda/PyPI 和 R 源文件，这套配置不是离线环境镜像。版本匹配不等同于逐文件审计安装内容。

## 维护

在这里修改配置和锁文件，再通过部署工具同步公共运行目录。更新补装版本前验证源文件及实际安装结果。修改环境文件或共同补装程序后，重新计算 `bundle.json` 的对应 SHA256，运行环境审计、选定分支依赖检查和 16 个 skill 的 fixture 流程；按仓库规则提交、推送。不要直接修改 `.pixi/` 文件来替代规范配置。
