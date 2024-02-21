# Documentation folder

The following documentation is available

| File                         | Description |
|------------------------------|-------------|
|CONTAINERS.md                 | How to install / use the SneakerNet container |
|INSTALL.advanced.md           | Some older installation notes; might have some esoteric details |
|INSTALL.md                    | Installation instructions |
|INSTALL.taxonProperties.md    | How to add a new organism so that SneakerNet recognizes it and analyzes it appropriately |
|PLUGINSDEV.md                 | How to develop your own plugin |
|PLUGINS.md                    | The plugins catalog |
|SneakerNetInput.md            | SneakerNet input description |
|SneakerNetOutput.md           | SneakerNet output description |

## Directory structure

```plaintext
SneakerNet
├── config                     Modified config files for your installation
├── config.bak                 Original config files from repo
├── db                         Database that is built during runtime
│   ├── bed
│   ├── fasta
│   ├── pubmlst_download
│   ├── staramr
│   └── wgMLST
├── docs                       Documentation on SneakerNet
│   ├── images
│   └── plugins
├── example                    Some examples on running SneakerNet
│   ├── images
│   ├── inbox
│   └── M00123-18-001-test
├── lib
│   └── perl5                  Perl library and any libraries imported during installation
├── man                        Perl library documentation imported during installation
│   ├── man1
│   └── man3
├── paper                      JOSS manuscript
├── scripts                    Main SneakerNet scripts
│   └── dockerAliases          Shell scripts to invoke certain docker containers
├── SneakerNet.plugins         SneakerNet plugins as described in PLUGINS.md
│   ├── deprecated
│   └── helper                 Any helper scripts for the plugins
└── t                          Unit tests
    ├── futureTests
    ├── M00123-18-001-test
    ├── M00123-18-002-test
    └── M00123-20-001-sarscov2
```
