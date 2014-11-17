
class Arithmetic:


    def length(self,N):
        return (1+N)

    def lengthstr(self,N):
        return "%.1f"%self.length(N)

    def clusters(self,N):
        x = N
        y = 1
        res = [(x,y),]
        while x > (y+1):
            x -= 1
            y += 1
            res.append((x,y))
        return res
