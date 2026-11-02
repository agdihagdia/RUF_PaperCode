using PyCall
export set_dateAxis,setDateTick,get_dateFormat,set_language
"""
myFmt=get_dateFormat("%m-%d")
plt.gca().xaxis.set_major_formatter(myFmt)
"""
function set_language(lan="en_US.utf8")    
    locale=pyimport("locale")
    locale.setlocale(locale.LC_ALL,lan)
end
function get_dateFormat(formatString)
    mdates = pyimport("matplotlib.dates")
    myFmt = mdates.DateFormatter(formatString)
    return myFmt
end
function set_dateAxis(axis_xory=plt.gca().xaxis)
    mdates = pyimport("matplotlib.dates")
    locator = mdates.AutoDateLocator()
    formatter = mdates.ConciseDateFormatter(locator)
    axis_xory.set_major_locator(locator)
    axis_xory.set_major_formatter(formatter)
end

function setDateTick(ax=plt.gca();minticks=3, maxticks=10)
    mdates=pyimport("matplotlib.dates")
    
    locator = mdates.AutoDateLocator(minticks=minticks, maxticks=maxticks)
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
end