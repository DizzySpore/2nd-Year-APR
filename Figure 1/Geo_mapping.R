# Map a country with one administrative region highlighted and selected
# locations marked. No input data file required.
#
# Edit the CONFIGURATION section below to change the country, the highlighted
# region, the marked locations and the output file.

# Install required packages if not already installed (uncomment if needed)
# install.packages(c("ggplot2", "sf", "rnaturalearth", "rnaturalearthdata"))
# install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")

# ==== CONFIGURATION ========================================================

# Optional custom R library location. Leave as NULL to use the default
# library paths, or set it e.g. via an environment variable:
#   lib_path <- Sys.getenv("R_CUSTOM_LIB", unset = NA); if (is.na(lib_path)) lib_path <- NULL
lib_path <- NULL

country          <- "Malawi"     # Natural Earth country name
highlight_region <- "Chikwawa"   # admin-1 unit to highlight (matched case-insensitively)
highlight_fill   <- "red"

# Locations to mark (decimal degrees, WGS84; approximate town centres).
# nudge_x / nudge_y / hjust control label placement to avoid overlaps.
places <- data.frame(
  name      = c("Lilongwe (capital)", "Blantyre", "Chikwawa"),
  longitude = c(33.7741, 35.0058, 34.8000),
  latitude  = c(-13.9626, -15.7861, -16.0350),
  type      = c("Capital", "City / town", "City / town"),
  nudge_x   = c(0.15, 0.15, -0.15),
  nudge_y   = c(0.00, 0.10, -0.10),
  hjust     = c(0, 0, 1)
)

# Marker shape per location type (21 = circle, 23 = diamond, 22 = square)
type_shapes <- c("Capital" = 23, "City / town" = 21)

map_padding <- 0.3   # degrees added around the country's bounding box
output_file <- paste0(tolower(country), "_", tolower(highlight_region), "_map.png")
plot_width  <- 6     # inches
plot_height <- 9

# ==== PACKAGES =============================================================

library(ggplot2)
library(sf)
for (pkg in c("rnaturalearth", "rnaturalearthdata", "rnaturalearthhires")) {
  library(pkg, character.only = TRUE, lib.loc = lib_path)
}

# ==== BOUNDARIES ===========================================================

country_sf <- ne_countries(scale = "large", country = country, returnclass = "sf")
regions    <- ne_states(country = country, returnclass = "sf")

if (nrow(country_sf) == 0) stop("Country '", country, "' not found in Natural Earth.")

highlight_sf <- regions[grepl(highlight_region, regions$name, ignore.case = TRUE), ]

if (nrow(highlight_sf) == 0) {
  print(sort(regions$name))
  stop("Region '", highlight_region, "' not found - check the names printed above.")
}

# Map extent derived from the country outline
bb <- st_bbox(country_sf)
xlim <- c(bb["xmin"] - map_padding, bb["xmax"] + map_padding)
ylim <- c(bb["ymin"] - map_padding, bb["ymax"] + map_padding)

# ==== PLOT =================================================================

p <- ggplot() +
  geom_sf(data = regions, fill = "grey90", color = "grey60", linewidth = 0.2) +
  geom_sf(data = highlight_sf, fill = highlight_fill, color = "black",
          linewidth = 0.4, alpha = 0.8) +
  geom_sf(data = country_sf, fill = NA, color = "black", linewidth = 0.5) +
  geom_point(data = places,
             aes(x = longitude, y = latitude, shape = type),
             size = 3, fill = "gold", color = "black", stroke = 0.8) +
  geom_text(data = places,
            aes(x = longitude + nudge_x, y = latitude + nudge_y,
                label = name, hjust = hjust),
            size = 3.5, fontface = "bold") +
  scale_shape_manual(values = type_shapes, name = NULL) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = country,
       subtitle = paste(highlight_region, "highlighted"),
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, face = "italic"),
        legend.position = "bottom")

print(p)

ggsave(output_file, plot = p, width = plot_width, height = plot_height, dpi = 300)
cat("Map saved as '", output_file, "'\n", sep = "")
