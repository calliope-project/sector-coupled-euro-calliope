# Rules to sync to and from Euler

EULER_URL = "euler.ethz.ch"
EULER_BASE_DIR = "~/Develop/money-land/"
EULER_BUILD_DIR = EULER_BASE_DIR + "build/"
LOCAL_EULER_RESULTS = "./build/euler"
SYNCIGNORE = ".syncignore"


rule send:
    message: "Send changes to Euler"
    shell:
        "rsync -avzh --progress --delete -r . --exclude-from={SYNCIGNORE} {EULER_URL}:{EULER_BASE_DIR}"


rule receive:
    message: "Receive build changes from Euler"
    shell:
        "rsync -avzh --progress --delete -r --exclude-from=.syncignore-build {EULER_URL}:{EULER_BUILD_DIR} {LOCAL_EULER_RESULTS}"


rule clean_euler:
    message: "Clean results downloaded from Euler"
    shell:
        """
        rm -r {LOCAL_EULER_RESULTS}/*
        """
