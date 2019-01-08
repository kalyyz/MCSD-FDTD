% The Bi-complex number class for 3-D FDTD simulations
% A bicomplex number: z1+z2j2
% In this class, z1 and z2 are defined as 3-D tensors for 3-D FDTD simulations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BI CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef bi 
   properties
      z1;
      z2;
   end
   methods
      function self = bi(z1,z2);
          self.z1=z1;
          self.z2=z2;
      end
      
      function  out = subsref(self,index) % Indexing
            if  strcmp('()' ,index.type)
                %b = index.subs(2);
                %c = size(b{1});
                %out = bi([],[]);
                %out.z1 = builtin('subsref' ,self.z1,index);
                %out.z2 = builtin('subsref' ,self.z2,index);
                %out = b
                %a = builtin('subsref' ,self.z1,index);
                %b = builtin('subsref' ,self.z2,index);
                %out = bi(a,b);
                x = index.subs{1};
                y = index.subs{2};
                z = index.subs{3};
                %out = self;
                self.z1 = self.z1(x,y,z);
                self.z2 = self.z2(x,y,z);
                out = self;
                
            elseif  strcmp('.' ,index.type)
                out = eval(['self.' ,index.subs]);
            end
          %{  
          elseif(length(index)==2)
            if  strcmp('()' ,index(1).type)
                if  strcmp(index(2).subs,'z1' )
                %out = bi([],[]);
                out = builtin('subsref' ,self.z1,index(1));
                
                elseif strcmp(index(2).subs,'z2' )
                out = builtin('subsref' ,self.z2,index(1));
                end
                
            end
          end
            %}
      end
      
      function  self = subsasgn(self,index,value) % Asigning
           %if(length(index)==1)
            if  strcmp('()' ,index.type)
                %self = bi([],[]);
                self.z1 = builtin('subsasgn' ,self.z1,index,value.z1);
                self.z2 = builtin('subsasgn' ,self.z2,index,value.z2);
                %x = index.subs{1};
                %y = index.subs{2};
                %z = index.subs{3};
                
                %self.z1(x,y,z) = value.z1;
                %self.z2(x,y,z) = value.z2;
            end
      end
   
      function out = plus(a,b);
          if (isobject(a)==1 && isobject(b)==1)
          out  = bi(a.z1+b.z1,a.z2+b.z2);    
          else
              if(isobject(b)==1)
                a = bi(a,0);    
                out  = bi(a.z1+b.z1,a.z2+b.z2);    
              else
                b = bi(b,0);    
                out  = bi(a.z1+b.z1,a.z2+b.z2);  
              end
          end
      end
      
      function out = minus(a,b);
          if (isobject(a)==1 && isobject(b)==1)
          out  = bi(a.z1-b.z1,a.z2-b.z2);    
          else
            if(isobject(b)==1)
                a = bi(a,0);    
                out  = bi(a.z1-b.z1,a.z2-b.z2);     
            else
                b = bi(b,0);
                out  = bi(a.z1-b.z1,a.z2-b.z2);     
            end
          end
      end
      
      function out = mtimes(a,b);
          if (isobject(a)==1 && isobject(b)==1)
            m = a.z1.*b.z1-a.z2.*b.z2;
            n = a.z1.*b.z2+a.z2.*b.z1;
            out = bi(m,n);
          else
              if(isobject(b)==1)
                a = bi(a,0)    
                b
                a.z1
                m = (a.z1).*(b.z1)-(a.z2).*(b.z2);
                n = a.z1.*b.z2+a.z2.*b.z1;
                out = bi(m,n);
              else
                b = bi(b,0);
                m = a.z1.*b.z1-a.z2.*b.z2;
                n = a.z1.*b.z2+a.z2.*b.z1;
                out = bi(m,n);
              end
          end
      end
      
      function out = mrdivide(a,b);
          if (isobject(a)==1&&isobject(b)==1)
          m = (a.z1.*b.z1+a.z2.*b.z2)./((b.z1).^2+(b.z2).^2);
          n = (a.z2.*b.z1-a.z1.*b.z2)./((b.z1).^2+(b.z2).^2);
          out = bi(m,n);
          else
              if(isobject(b)==1)
                  a = bi(a,0);
                  m = (a.z1.*b.z1+a.z2.*b.z2)./((b.z1).^2+(b.z2).^2);
                  n = (a.z2.*b.z1-a.z1.*b.z2)./((b.z1).^2+(b.z2).^2);
                  out = bi(m,n);
              else
                  b = bi(b,0);
                  m = (a.z1.*b.z1+a.z2.*b.z2)./((b.z1).^2+(b.z2).^2);
                  n = (a.z2.*b.z1-a.z1.*b.z2)./((b.z1).^2+(b.z2).^2);
                  out = bi(m,n);
              end
          end
                  
      end
      
      function out = sec(obj,x1,x2,y1,y2)
          m = obj.z1(x1:x2,y1:y2); 
          n = obj.z2(x1:x2,y1:y2);
          out = bi(m,n);
      end
      
       function update(obj,obj1,x1,x2,y1,y2)
          obj.z1(x1:x2,y1:y2) = obj1.z1(x1:x2,y1:y2); 
          obj.z2(x1:x2,y1:y2) = obj1.z2(x1:x2,y1:y2);
       end
      
       function assign(obj,x,y,obj1,x1,y1);
           obj.z1(x,y) = obj1.z1(x1,y1);
           obj.z2(x,y) = obj1.z2(x1,y1);
       end
       
       function eq(obj,obj1)
           obj.z1 = obj1.z1;
           obj.z2 = obj1.z2;
       end
       
       function out = u(obj,x,y)
            out = bi([],[]);
            out.z1 = obj.z1(x,y);
            out.z2 = obj.z2(x,y);
       end
      
       function out = real(obj)
           out = real(obj.z1);
       end
       
        function out = im1(obj)
           out = imag(obj.z1);
        end
       
         function out = im2(obj)
           out = real(obj.z2);
         end
     
        function out = im12(obj)
           out = imag(obj.z2);
        end
       
        function out = sqrtbi(obj);
            out = bi([],[]);
            out.z1 = sqrt( ( obj.z1 + sqrt( (obj.z1).^2 + (obj.z2).^2 ))/2 );
            out.z2 = sqrt( ( -obj.z1 + sqrt( (obj.z1).^2 + (obj.z2).^2 ))/2 );
        end
        
        function out = power(self,other) % Elementwise power
             sizes = size(self);
                z_1 = zeros(sizes);
                z_2 = zeros(sizes);

                for i = 1:length(z_1(:))
                 sr.type = '()';
                 sr.subs = {i};
                 r = 1;
                 theta = argc(subsref(self,sr));
                 z_1(i) = r^other*cos(other*theta);
                 z_2(i) = r^other*sin(other*theta);
                end
                out = bi([],[]);
                out.z1 = z_1;
                out.z2 = z_2;

        end
          
        function out = sin(self) % sin
            out = bi([],[]);
            out.z1=cosh(self.z2).*sin(self.z1);
            out.z2=sinh(self.z2).*cos(self.z1);
        end

        function out = cos(self) % cos
             out = bi([],[]);
             out.z1=cosh(self.z2).*cos(self.z1);
             out.z2=-sinh(self.z2).*sin(self.z1);
        end
        
        

       
   end
end
      
