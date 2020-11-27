classdef TABLE
  properties
    index
    data
  end
  methods
    %-----------------------------------------------------------constructor
    function this=TABLE(index,data)
      this.index = 0;
      this.data  = zeros(1,2);
      if nargin == 1
        this.index = index;
      elseif nargin == 2
        this.index = index;
        this.data  = data;
      end
    end
    %----------------------------------------------------------------------
    function this=set.index(this,index)
      this.index = index;
    end
    %----------------------------------------------------------------------
    function index=get.index(this)
      index =this.Index;
    end
    %----------------------------------------------------------------------
    function this=set.data(this,data)
      this.data = data;
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    function [val]=getTableVal(this,x)
      val    = 0.0;
      dimTab = size(this.data,1);
      if dimTab==0.0
          return
      end
      if x<=this.data(1,1)
          val = this.data(1,2);
          return
      end
      if x>=this.data(dimTab,1)
          val = this.data(dimTab,2);
          return
      end

      for i=2:dimTab
          if (x>this.data(i-1,1))&&(x<=this.data(i,1))
              x1  = this.data(i-1,1);
              x2  = this.data(i  ,1);
              y1  = this.data(i-1,2);
              y2  = this.data(i  ,2);
              val = y2-(y2-y1)*(x2-x)/(x2-x1);
              break;
          end
      end
    end
    %----------------------------------------------------------------------
    function [val] = getTableDVal(this, x)
      val    = 0.0;
      dimTab = size(this.data,1);
      if dimTab==0.0
        return
      end
      if x <= this.data(1,1)
        val = 0.0;
        return
      end
      if x >= this.data(dimTab, 1)
        if dimTab == 1
          val = 0.0;
          return
        end
        x2 = this.data(dimTab    , 1);
        x1 = this.data(dimTab - 1, 1);
        y2 = this.data(dimTab    , 2);
        y1 = this.data(dimTab - 1, 2);
        val = (y2 - y1) / (x2 - x1);
        return
      end

      for i=2:dimTab
          if (x>this.data(i-1,1))&&(x<=this.data(i,1))
              x1  = this.data(i-1,1);
              x2  = this.data(i  ,1);
              y1  = this.data(i-1,2);
              y2  = this.data(i  ,2);
              val = (y2 - y1) / (x2 - x1);
              break;
          end
      end
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
  end
end