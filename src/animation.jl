using Javis
using Interpolations

#object constructors
function ground(args...)
    background("white")
    sethue("black")
end
function projectile(xMin, yMin, length, radius)
    sethue("red")
    rect(xMin, yMin, length, 2radius, :fill)
end
function drawtext(xMin, yMin, string, size)
    fontsize(size)
    Javis.text(string,[xMin,yMin])
end
function barrel(xMin, length, bore)
    sethue("grey")
    rect(xMin, bore/2, length, 5, :fill)
    rect(xMin, -bore/2-5, length, 5, :fill)
end
function coil(position, length, radius)
    setcolor(sethue("grey")..., 0.5)
    rect(position-length, -radius-20, length, 2*radius+40, :fill)
end
function magcoil(position, length, radius)
    setcolor(sethue("red")..., 0.5)
    rect(position-length, -radius-20, length, 2*radius+40, :fill)
end
function revmagcoil(position, length, radius)
    setcolor(sethue("blue")..., 0.5)
    rect(position-length, -radius-20, length, 2*radius+40, :fill)
end
function coilevents(i)
    coilMid = ustrip(scenario.coils[i].location + scenario.coils[i].length/2 - scenario.proj.physical.length*1.5)
    coilMidTime = LinearInterpolation(sln[1,:], sln.t).(coilMid)
    coilTimes=[ustrip(scenario.eventTimes.entersActiveZone[i]),coilMidTime,ustrip(scenario.eventTimes.exitsActiveZone[i])]
    return coilTimes
end

#animation parameters
vidlength = 1000
vidheight = 500
video = Video(vidlength, vidheight)
zero = -vidlength/2
framecount = 250
trimtowindow = true

#scale and trim data to animation
animscale=vidlength/(scenario.barrel.length+2*scenario.proj.physical.length)
rangemax=length(sln.t)
if trimtowindow
    rangemax = 1
    while ustrip(sln[1, rangemax]*animscale) < vidlength
        global rangemax+=1
    end
end
posfr=ustrip(LinearInterpolation(sln.t, sln[1,:]).(range(0,sln.t[rangemax],length=framecount))*animscale)
timescale=framecount/sln.t[rangemax]

#generate objects
anim_background = Background(1:framecount, ground)
Object((args...) -> barrel(zero, scenario.barrel.length*animscale, 50))
for n in 1:length(scenario.coils)
    coilFrames = [Int(round(timescale * i)) for i in coilevents(n)]
    currentcoil=scenario.coils[n]
    Object(1:coilFrames[1], (args...) -> coil(currentcoil.location*animscale+zero, currentcoil.length*animscale-5, 25))
    Object((coilFrames[1]==0 ? 1 : coilFrames[1]+1):coilFrames[2], (args...) -> 
        magcoil(currentcoil.location*animscale+zero, currentcoil.length*animscale-5, 25))
    Object((coilFrames[2]==0 ? 1 : coilFrames[2]+1):coilFrames[3], (args...) -> 
        revmagcoil(currentcoil.location*animscale+zero, currentcoil.length*animscale-5, 25))
    Object((coilFrames[3]+1):framecount, (args...) -> 
        coil(currentcoil.location*animscale+zero, currentcoil.length*animscale-5, 25))
end
proj = Object(1:framecount, (args...) -> projectile(zero,-25, scenario.proj.physical.length*animscale, 25))
for i in 2:length(posfr)
    act!(proj,Action(i:i, anim_translate(posfr[i]-posfr[i-1], 0),),)
    Object(i:i, (args...) -> Javis.text(string(i); halign = :center))
end
render(video, pathname="animation.gif", framerate = 25)