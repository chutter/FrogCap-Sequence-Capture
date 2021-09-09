# Package
version       = "0.1.0"
author        = "Kerry Cobb & Carl Hutter"
description   = "FrogCap Documentation"
license       = "MIT"
skipDirs      = @["book"]
bin           = @["nbook"]
binDir        = "bin"

# Dependencies
requires "nim >= 1.2.0"
requires "https://github.com/pietroppeter/nimibook"


# task genbook, "build book":
#   exec("nimble build -d:release")
#   exec("./bin/nbook init")
#   exec("./bin/nbook build")