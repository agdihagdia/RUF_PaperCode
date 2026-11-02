export addMultiText
function addMultiText(ha,x_all,y_all,str_all;hori_ali=:right,fontsize=8,color_all=["k" for i in 1:length(str_all)])
    for (n,str) in enumerate(str_all)
        ha.text(x_all[n],y_all[n],str,transform=ha.transAxes,verticalalignment="center", horizontalalignment=hori_ali,fontsize=fontsize,color=color_all[n])
    end
end