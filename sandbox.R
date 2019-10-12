#A sandbox. The extension helps my text editor (Visual Studio Code) to format the code.
txt <- ""
count_down <- 3
position <- 1
while (count_down >= 0) {
    if (count_down > 0) {
        txt[position] <- as.character(count_down)
        #print(as.character(count_down))
    } else {
        txt[position] <- "Lift Off!"
    }
    count_down <- count_down - 1
    position <- position + 1
}

print(txt)

countDown <- function(x=3) {
    txt <- ""
    start <- x
    position <- 1
    while (start >= 0) {
        if (start > 0) {
            txt[position] <- as.character(start)
        } else {
            txt[position] <- "Lift Off!"
        }
        start <- start - 1
        position <- position + 1
    }
    return(txt)
}

myLifeDays <- function(birthday, lday) {
    if (missing(birthday)) {
        print ("Enter your birthday as a string in \"YYYY-MM-DD\" format.")
        return()
    }
    bd <- strptime(birthday, "%Y-%m-%d") # convert string to time
    now <- format(Sys.time(), "%Y-%m-%d") # convert "now" to time
    diff <- round(as.numeric(difftime(now, bd, unit="days")))
    print(sprintf("This date was %d days ago.", diff))

    futureDay <- format(bd + as.difftime(lday, unit="days"), "%Y-%m-%d") #Consulted ?difftime
    print(sprintf("Celebrate your life day on: %s", futureDay))
    return()
}