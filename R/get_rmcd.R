get_rmcd <- function(data, rmcd = "dat/ODP208_1267_rmcd.csv") {
  rmcd <- readr::read_csv(rmcd) |>
    separate(label, into = c("sitehole", "coretype", "Sec"),
               sep = "-", remove = FALSE) |>
    separate(sitehole, into = c("Site", "H"), sep = -1) |>
    separate(coretype, into = c("Core", "T"), sep = -1) |>
    # we do not rename the interval, may not be the same as in the data!
    # rename the CC sections into 7, the naming convention in the MS data
    mutate(Sec = ifelse(Sec == "7", "7", Sec),
           Sec = ifelse(Sec == "cc", "C", Sec)) |>
    mutate(diff = depth_rmcd - depth_mbsf, .after = depth_rmcd) |>
    mutate(diff2 = depth_rmcd2 - depth_mbsf2, .after = depth_rmcd2) |>
    mutate(row = 1:n())

  # the right side of the splice table only
  rmcd2 <- rmcd |>
    select(label, link, label2, interval2, depth_mbsf2, depth_rmcd2, diff2, row) |>
    separate(label2, into = c("sitehole", "coretype", "Sec"),
               sep = "-", remove = FALSE) |>
    separate(sitehole, into = c("Site", "H"), sep = -1) |>
    separate(coretype, into = c("Core", "T"), sep = -1) |>
    # we do not rename the interval, may not be the same as in the data!
    mutate(Sec = ifelse(Sec == "7", "7", Sec),
           Sec = ifelse(Sec == "cc", "C", Sec))

  out <- data |>
    tidylog::left_join(rmcd |>
                     # make the types the same
                     mutate(across(c(Site, Core), parse_double)) |>
                     # do NOT match by section, only by core!
                     rename(section = Sec) |>
                     select(top = label, to = label2,
                            Site, H, Core, T, section, interval,
                            depth_mbsf, depth_rmcd, diff, row)) |>
    # add the right-hand side of the splice table
    tidylog::left_join(rmcd2 |>
                       mutate(across(c(Site, Core), parse_double)) |>
                       rename(section2 = Sec) |>
                       select(from = label, bot = label2,
                              Site, H, Core, T, section2, interval2,
                              depth_mbsf2, depth_rmcd2, diff2, row2 = row)) |>
  mutate(my_rmcd = case_when(
  (Sec <= section) | ((Sec == section) & (`Top (cm)` <= interval)) ~
    `Depth (mbsf)` + diff,
  (Sec >= section2) | ((Sec == section2) & (`Top (cm)` >= interval2)) ~
    `Depth (mbsf)` + diff2,
  TRUE ~ NA_real_)) |>
    mutate(on_splice = (Sec < section | ((Sec == section) &
                                         (`Top (cm)` <= interval))) &
             (Sec > section2 | ((Sec == section2) &
                                `Top (cm)` >= interval2))) |>
    mutate(on_splice = ifelse(is.na(on_splice), FALSE, on_splice))

  return(out)
}
