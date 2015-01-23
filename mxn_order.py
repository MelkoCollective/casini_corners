class Arithmetic:


    def length(self,N):
        return (1+N)

    def lengthstr(self,N):
        return "%.1f"%self.length(N)

    def clusters(self,N):
        x = N-1
        y = 1
        res = [(x,y),]
        while x > (y+1):
            x -= 1
            y += 1
            res.append((x,y))
        return res

class Geometric:

    def clusters(self,N):
        res = [(N,1),]
        i = 2
        sqrtN = N ** (0.5)
        while i < sqrtN:
            if N % i == 0:
                res.append((N/i,i))
            i += 1
        if i*i == N:
            res.append((int(sqrtN),int(sqrtN)))

        return res

