function z=calcarea(w,h,corx,cory,uppx,uppy,ddx,ddy)
    arect = (w-corx)*(h-cory);
    atu = (h-cory)*(corx-uppx)/2;
    atd = (w-corx)*(cory-ddy)/2;
    z = arect+atu+atd;
end
