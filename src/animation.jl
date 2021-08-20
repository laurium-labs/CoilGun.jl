using Javis
using Interpolations

#object constructors
function ground(args...)
    background("white")
    sethue("black")
end
function projectile(xMin, yMin, length, radius)
    rbBlend = blend(Point(xMin,0), Point(xMin+length,0), "red", "blue")
    setblend(rbBlend)
    rect(xMin, yMin, length, 2radius, :fill)
end
function drawtext(xMin, yMin, string, size)
    fontsize(size)
    Javis.text(string,[xMin,yMin])
end
function barrel(xMin, length, bore)
    sethue("grey")
    rect(xMin, bore / 2, length, 5, :fill)
    rect(xMin, -bore / 2 - 5, length, 5, :fill)
end
function coil(position, length, radius)
    setcolor(sethue("grey")..., 0.5)
    rect(position - length, -radius - 20, length, 2 * radius + 40, :fill)
end
function magcoil(position, length, radius)
    xMin = position - length
    rbBlend = blend(Point(xMin, 0), Point(xMin + length, 0))
    addstop(rbBlend, 0, setcolor(sethue("red")..., 0.7))
    addstop(rbBlend, 1, setcolor(sethue("blue")..., 0.7))
    setblend(rbBlend)
    rect(xMin, -radius - 20, length, 2 * radius + 40, :fill)
end
function revmagcoil(position, length, radius)
    xMin = position - length
    brBlend = blend(Point(xMin, 0), Point(xMin + length, 0))
    addstop(brBlend, 0, setcolor(sethue("blue")..., 0.7))
    addstop(brBlend, 1, setcolor(sethue("red")..., 0.7))
    setblend(brBlend)
    rect(xMin, -radius - 20, length, 2 * radius + 40, :fill)
end
function magneticFieldLine(ellipseCenter, seperationDist, height, hueIntensity)
    setcolor(sethue("green1")..., hueIntensity)
    ellipse(
        Point(ellipseCenter - seperationDist/2, height), 
        Point(ellipseCenter + seperationDist/2, height), 
        Point(ellipseCenter, 0), 
        :stroke
    )
end
function magneticFieldLines(coilCenter, coilLength, coilRadius)
    c = coilCenter .- coilLength/2
    rad  = coilRadius + 20
    adj  = 0.75
    for scaler in [-adj,adj]
        magneticFieldLine(c, coilLength, scaler * rad * 1.00, 0.6)
        magneticFieldLine(c, coilLength, scaler * rad * 1.25, 0.4)
        magneticFieldLine(c, coilLength, scaler * rad * 1.50, 0.3)
        magneticFieldLine(c, coilLength, scaler * rad * 1.75, 0.1)
    end
end
function coilevents(i)
    coilMid = ustrip(scenario.coils[i].location - scenario.coils[i].length/2)# + scenario.proj.physical.length * 1.5)
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
animscale = vidlength / (scenario.barrel.length + 2 * scenario.proj.physical.length)
rangemax = if trimtowindow
    findlast(pos -> ustrip(pos * animscale) < vidlength, sln[1,:]) + 1
else
    length(sln.t)
end
posfr = ustrip(LinearInterpolation(sln.t, sln[1,:]).(range(0,sln.t[rangemax],length = framecount)) * animscale)
timescale = framecount / sln.t[rangemax]

#generate objects
anim_background = Background(1:framecount, ground)
Object((args...) -> barrel(zero, (scenario.barrel.length + scenario.coils[1].length) * animscale, 50))
for n in eachindex(scenario.coils)
    coilFrames = [round(Int, timescale * t) for t in coilevents(n)]
    @show coilFrames
    currentcoil = scenario.coils[n]

    coilPosition = currentcoil.location * animscale + (zero + currentcoil.length*animscale)
    coilLength = currentcoil.length * animscale - 5
    coilRadius = 25
    inputArgs = (coilPosition, coilLength, coilRadius)
    currentIn = (coilFrames[1] + 1):coilFrames[3]
    currentOut = (coilFrames[3] + 1):(n <= lastindex(scenario.coils) - 2 ? round(Int, timescale * coilevents(n + 2)[end]) : framecount)

    Object(currentIn, (args...) -> magneticFieldLines(inputArgs...))
    Object(currentOut, (args...) -> magneticFieldLines(inputArgs...))
    Object(1:framecount, (args...) -> coil(inputArgs...))
    Object(currentIn, (args...) -> magcoil(inputArgs...))
    Object(currentOut, (args...) -> revmagcoil(inputArgs...))
end

proj = Object(1:framecount, (args...) -> projectile(zero, -25, scenario.proj.physical.length * animscale, 25))

prevPoint = O
for (idx, pos) in enumerate(posfr)
    act!(proj, Action(idx:idx, anim_translate(prevPoint, Point(pos,0))))
    # Object(idx:idx, (args...) -> Javis.text(string(idx); halign = :center))
    global prevPoint = Point(pos,0)
end
render(video, pathname = "animation.gif", framerate = 25)