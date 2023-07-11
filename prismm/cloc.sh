find . -name '*.py' -or -name '*.sh' -or -name '*.sage' -or -name '*.R' | xargs wc -l
find . -name '*.py' -or -name '*.sh' -or -name '*.sage' -or -name '*.R' -exec cat {} \; | wc -w


