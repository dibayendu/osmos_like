# osmos_like

This is a game which is like osmos. I am writing this code to practise OpenGL using C and C++

- This is game which is like OSMOS.
- This is just a way to make a game like OSMOS and not to replicate it.
- This project should not be commercially
- The code is writtern in C and C++
- using OpenGL framework and libraries.


Assumptions/Notes:
----

1. Conserved momentum of the circular arena.
2. In uniform grid if cell size is smalled than mote diameter, I am using brute force in that scenario.
3. Starting the game with 100 motes.
4. In my code, I have called mote as particle too where both means the same.
5. The player is the one which has the texture of the lookalike texture of osmos player mote.
6. The player is always at position 0,0 of the screen.
7. Motes other than the player mote are given random colours.
8. Code indentation is done keeping in view it will be read in a full screen vim environment.
9. Many printf statements are commented out which was used in debugging and can be reused for testing 
    e.g. checking if the momentum of the environment is conserved.
10. Once game is over there will not be any mouse interactions.
11. If the player mote is absorbed, the camera will be focused at the center of the arena.
 

Version
----
1.1


Feature implemented/Not implemented/bug:
----
1. Shaders are not implemented.
2. If any mote other than the player mote, has more than or equal to half of the mass of the 
   environment, then the game cannot be won, but this message is not shown to the player,
   keeping the player still alive and not ending the game, even though it cannot be won.


COMMANDS:
----
1. Compile  : make
2. Run      : "./assignment_3_application 1"    # 1 = square; 2 = circle
            : make square       # starts a square arena
            : make circle       # starts a circular arena
3. Clean    : make clean
4. Archive  : make archive


FUNCTIONALITY:
----
keyboard
----
    1.  Esc/q   : quit
    2.  w       : wireframe/solid view
    3.  d       : bruteForce/uniformGrid collision detection
    4.  s       : start/pause the game
    5.  u       : restart the game

mouse
----
    1. leftclick    : release motes in a particular direction and the player mote gains velocity.
    2. rightclick   : zooms in/out
    

OTHERS:
----
1. course website tutorial examples
3. Osmos game textures
4. google images textures
