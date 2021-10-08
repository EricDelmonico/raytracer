using Microsoft.Xna.Framework;
using Microsoft.Xna.Framework.Graphics;
using Microsoft.Xna.Framework.Input;
using System.Threading;

namespace raytracer
{
    /// <summary>
    /// This is the main type for your game.
    /// </summary>
    public class Game1 : Game
    {
        GraphicsDeviceManager graphics;
        SpriteBatch spriteBatch;

        // the texture that will render the 
        // raytraced scene to the window
        private Texture2D render;

        private float FPS = 0;
        private SpriteFont arial16;
        private Point screenSz;

        private bool UI = true;

        private KeyboardState kbState;
        private KeyboardState prev_kbState;

        public Game1()
        {
            graphics = new GraphicsDeviceManager(this);
            float scale = 2;
            screenSz = new Point((int)(Raytracing.IMG_WIDTH * scale), (int)(Raytracing.IMG_HEIGHT * scale));
            graphics.PreferredBackBufferWidth = screenSz.X;
            graphics.PreferredBackBufferHeight = screenSz.Y;
            Content.RootDirectory = "Content";
        }

        /// <summary>
        /// Allows the game to perform any initialization it needs to before starting to run.
        /// This is where it can query for any required services and load any non-graphic
        /// related content.  Calling base.Initialize will enumerate through any components
        /// and initialize them as well.
        /// </summary>
        protected override void Initialize()
        {
            // TODO: Add your initialization logic here
            render = new Texture2D(GraphicsDevice, 
                                   Raytracing.IMG_WIDTH, 
                                   Raytracing.IMG_HEIGHT);

            Raytracing.Init();

            // get the texture's data from Raytracing class
            ThreadPool.QueueUserWorkItem(new WaitCallback((object state) => Raytracing.GetTextureDataThreadedByScanline()));

            base.Initialize();
        }

        /// <summary>
        /// LoadContent will be called once per game and is the place to load
        /// all of your content.
        /// </summary>
        protected override void LoadContent()
        {
            // Create a new SpriteBatch, which can be used to draw textures.
            spriteBatch = new SpriteBatch(GraphicsDevice);

            // TODO: use this.Content to load your game content here
            arial16 = Content.Load<SpriteFont>("arial16");
        }

        /// <summary>
        /// UnloadContent will be called once per game and is the place to unload
        /// game-specific content.
        /// </summary>
        protected override void UnloadContent()
        {
            // TODO: Unload any non ContentManager content here
        }

        /// <summary>
        /// Allows the game to run logic such as updating the world,
        /// checking for collisions, gathering input, and playing audio.
        /// </summary>
        /// <param name="gameTime">Provides a snapshot of timing values.</param>
        protected override void Update(GameTime gameTime)
        {
            if (GamePad.GetState(PlayerIndex.One).Buttons.Back == ButtonState.Pressed || Keyboard.GetState().IsKeyDown(Keys.Escape))
                Exit();

            kbState = Keyboard.GetState();
            if (kbState.IsKeyDown(Keys.U) && !prev_kbState.IsKeyDown(Keys.U))
            {
                UI = !UI;
            }

            render.SetData(Raytracing.data);

            prev_kbState = kbState;
            // TODO: Add your update logic here
            base.Update(gameTime);
        }

        /// <summary>
        /// This is called when the game should draw itself.
        /// </summary>
        /// <param name="gameTime">Provides a snapshot of timing values.</param>
        protected override void Draw(GameTime gameTime)
        {
            GraphicsDevice.Clear(Color.Black);

            // TODO: Add your drawing code here
            spriteBatch.Begin(SpriteSortMode.Deferred, null, SamplerState.PointClamp);

            spriteBatch.Draw(render, new Rectangle(Point.Zero, screenSz), Color.White);

            // draw framerate at the top of the screen
            FPS = (float)(1.0 / gameTime.ElapsedGameTime.TotalSeconds);

            if (UI)
            {
                spriteBatch.DrawString(arial16, $"FPS: {FPS}", Vector2.Zero, Color.Red);
                spriteBatch.DrawString(arial16, $"Scanlines Remaining: {Raytracing.scanlinesRemaining}", new Vector2(0, 21), Color.Red);
                spriteBatch.DrawString(arial16, $"Free Threads (0-4): {Raytracing.maxThreads}", new Vector2(0, 42), Color.Red);
            }

            spriteBatch.End();
            base.Draw(gameTime);
        }
    }
}
