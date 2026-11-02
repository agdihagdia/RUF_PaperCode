module FigureStyle

export FIGURE_WIDTH,
       PANEL_TITLE_SIZE,
       AXIS_LABEL_SIZE,
       TICK_SIZE,
       ANNOTATION_SIZE,
       PANEL_LEGEND_SIZE,
       FIGURE_LEGEND_SIZE,
       ROW_LABEL_SIZE,
       AXIS_LINEWIDTH,
       DATA_LINEWIDTH,
       GRID_LINEWIDTH,
       OBS_COLOR,
       AI_COLOR,
       AURORA_COLOR,
       RUF_COLOR,
       NEUTRAL_COLOR,
       font_properties,
       configure_matplotlib!

const FIGURE_WIDTH = 6.85

# Use a minimum figure text size of 8 pt for readability.
# The hierarchy below keeps dense multi-panel figures readable without letting
# labels dominate the plotting area.
const PANEL_TITLE_SIZE = 9.5
const AXIS_LABEL_SIZE = 9.0
const TICK_SIZE = 8.0
const ANNOTATION_SIZE = 8.0
const PANEL_LEGEND_SIZE = 8.0
const FIGURE_LEGEND_SIZE = 8.5
const ROW_LABEL_SIZE = 9.5

const AXIS_LINEWIDTH = 0.7
const DATA_LINEWIDTH = 1.1
const GRID_LINEWIDTH = 0.5

# Preserve the established paper palette while pairing colors with line styles
# and markers so interpretation does not rely on color alone.
const OBS_COLOR = "deepskyblue"
const AI_COLOR = "#2ca02c"
const AURORA_COLOR = "#F58518"
const RUF_COLOR = "#C44E52"
const NEUTRAL_COLOR = "deepskyblue"

function font_properties(font_manager, regular_path, bold_path)
    return (
        panel_title = font_manager.FontProperties(fname=bold_path, size=PANEL_TITLE_SIZE),
        axis_label = font_manager.FontProperties(fname=regular_path, size=AXIS_LABEL_SIZE),
        tick = font_manager.FontProperties(fname=regular_path, size=TICK_SIZE),
        annotation = font_manager.FontProperties(fname=regular_path, size=ANNOTATION_SIZE),
        panel_legend = font_manager.FontProperties(fname=regular_path, size=PANEL_LEGEND_SIZE),
        figure_legend = font_manager.FontProperties(fname=regular_path, size=FIGURE_LEGEND_SIZE),
        row_label = font_manager.FontProperties(fname=bold_path, size=ROW_LABEL_SIZE),
    )
end

function configure_matplotlib!(matplotlib)
    # PyCall converts rcParams to a Julia Dict when accessed as a property, so
    # mutating that converted value does not change Matplotlib.  Use the public
    # rc() interface to update the live Python configuration instead.
    matplotlib.rc("font", family="DejaVu Sans", size=TICK_SIZE)
    matplotlib.rc("axes", titlesize=PANEL_TITLE_SIZE,
        labelsize=AXIS_LABEL_SIZE, linewidth=AXIS_LINEWIDTH)
    matplotlib.rc("xtick", labelsize=TICK_SIZE)
    matplotlib.rc("ytick", labelsize=TICK_SIZE)
    matplotlib.rc("legend", fontsize=PANEL_LEGEND_SIZE)
    matplotlib.rc("lines", linewidth=DATA_LINEWIDTH)
    matplotlib.rc("pdf", fonttype=42)
    matplotlib.rc("ps", fonttype=42)
    matplotlib.rc("savefig", dpi=600)
    matplotlib.rc("figure", dpi=600)
    return nothing
end

end
