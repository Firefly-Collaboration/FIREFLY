document.addEventListener('DOMContentLoaded', function() {
  console.log('Page loaded');

  const topbar = document.getElementById('topbar');
  const menuIcon = document.getElementById('menuIcon');
  const menuDrawer = document.getElementById('menuDrawer');
  const githubButton = document.getElementById('githubButton');

  // Menu drawer toggle
  if (menuIcon && menuDrawer) {
    menuIcon.addEventListener('click', function() {
      menuDrawer.classList.toggle('open');
    });

    // Optional: close menu when clicking outside
    document.addEventListener('click', function(e) {
      if (!menuDrawer.contains(e.target) && !menuIcon.contains(e.target)) {
        menuDrawer.classList.remove('open');
      }
    });
  }

  // GitHub redirect
  if (githubButton) {
    githubButton.addEventListener('click', function() {
      window.open('https://github.com/FireflySpectra/firefly_release', '_blank');
    });
  }

  // Python version guide redirect (if present on page)
  const pythonButton = document.querySelector('.btn-python');
  if (pythonButton) {
    pythonButton.addEventListener('click', function() {
      window.open('https://devguide.python.org/versions/', '_blank');
    });
  }

  // Hero parallax background and depth lift with subtle zoom
  const heroShell = document.querySelector('.ffx-hero-shell');
  const heroLayer = document.querySelector('.ffx-hero-layer');
  if (heroShell) {
    const heroHeight = heroShell.offsetHeight || 1;
    let ticking = false;
    let lastScroll = 0;

    const updateParallax = () => {
      const scrollY = window.scrollY || window.pageYOffset || 0;
      
      // Throttle updates if scroll hasn't changed much
      if (Math.abs(scrollY - lastScroll) < 5) {
        ticking = false;
        return;
      }
      lastScroll = scrollY;
      
      const clamped = Math.min(scrollY, heroHeight * 1.2);
      const bgShift = clamped * 0.08;
      const zoomAmount = 1 + (clamped / heroHeight) * 0.03;
      
      heroShell.style.setProperty('--ffx-hero-shift', `${bgShift}px`);
      heroShell.style.setProperty('--ffx-hero-zoom', zoomAmount);
      
      if (heroLayer) {
        heroLayer.style.setProperty('--ffx-hero-content-shift', `${bgShift * 0.25}px`);
      }
      ticking = false;
    };

    const onScroll = () => {
      if (!ticking) {
        window.requestAnimationFrame(updateParallax);
        ticking = true;
      }
    };

    updateParallax();
    window.addEventListener('scroll', onScroll, { passive: true });
  }

  // Reveal on scroll
  const revealables = document.querySelectorAll('.ffx-reveal');
  if (revealables.length) {
    const io = new IntersectionObserver((entries) => {
      entries.forEach((entry) => {
        if (entry.isIntersecting) {
          entry.target.classList.add('ffx-visible');
        } else {
          entry.target.classList.remove('ffx-visible');
        }
      });
    }, { threshold: 0.1, rootMargin: '0px 0px -5% 0px' });

    revealables.forEach((el) => io.observe(el));
  }

  // Space sections with star generation and parallax
  const spaceSections = document.querySelectorAll('.ffx-space-section');
  if (spaceSections.length) {
    let ticking = false;
    let mouseX = 0.5;
    let mouseY = 0.5;

    // Generate stars for each section
    spaceSections.forEach((section) => {
      const backLayer = section.querySelector('.ffx-space-stars-back');
      const midLayer = section.querySelector('.ffx-space-stars-mid');
      
      // Generate 200 background stars
      if (backLayer) {
        for (let i = 0; i < 200; i++) {
          const star = document.createElement('div');
          star.className = 'ffx-space-star';
          const size = Math.random() * 1.5 + 0.5;
          star.style.width = `${size}px`;
          star.style.height = `${size}px`;
          star.style.top = `${Math.random() * 100}%`;
          star.style.left = `${Math.random() * 100}%`;
          star.style.opacity = Math.random() * 0.5 + 0.3;
          star.style.animationDelay = `${Math.random() * 3}s`;
          backLayer.appendChild(star);
        }
      }
      
      // Generate 100 foreground stars
      if (midLayer) {
        for (let i = 0; i < 100; i++) {
          const star = document.createElement('div');
          star.className = 'ffx-space-star';
          const size = Math.random() * 2 + 1;
          star.style.width = `${size}px`;
          star.style.height = `${size}px`;
          star.style.top = `${Math.random() * 100}%`;
          star.style.left = `${Math.random() * 100}%`;
          star.style.opacity = Math.random() * 0.6 + 0.4;
          star.style.animationDelay = `${Math.random() * 3}s`;
          star.style.animationDuration = `${Math.random() * 2 + 2}s`;
          midLayer.appendChild(star);
        }
      }
    });

    const updateSpace = () => {
      const vh = window.innerHeight || 1;
      const viewportCenter = vh / 2;
      
      spaceSections.forEach((section) => {
        const back = section.querySelector('.ffx-space-stars-back');
        const mid = section.querySelector('.ffx-space-stars-mid');
        
        // Update each content element individually
        const contentElements = section.querySelectorAll('.ffx-space-card, .ffx-space-text > *, .ffx-space-imagewrap, .ffx-models-head, .ffx-models-grid, .ffx-models-bottom, .table-container');
        contentElements.forEach((element) => {
          const rect = element.getBoundingClientRect();
          const elementCenter = rect.top + rect.height / 2;
          const distance = viewportCenter - elementCenter;
          const progress = Math.min(Math.max(1 - Math.abs(distance) / (vh * 0.5), 0), 1);
          element.style.setProperty('--ffx-motion', progress.toFixed(3));
        });

        // Subtle scroll-driven drift for stars based on section center
        const sectionRect = section.getBoundingClientRect();
        const sectionCenter = sectionRect.top + sectionRect.height / 2;
        const sectionDistance = viewportCenter - sectionCenter;
        const scrollDrift = Math.max(-1, Math.min(1, sectionDistance / (vh * 0.8)));
        
        // Reduced parallax intensity for smoother performance
        const parallaxX = (mouseX - 0.5) * 50;
        const parallaxY = (mouseY - 0.5) * 25;

        if (back) {
          back.style.transform = `translate3d(${parallaxX * 0.2}px, ${parallaxY * 0.2 + scrollDrift * 10}px, 0)`;
        }
        if (mid) {
          mid.style.transform = `translate3d(${parallaxX * 0.5}px, ${parallaxY * 0.5 + scrollDrift * 18}px, 0)`;
        }
      });
      ticking = false;
    };

    const onScroll = () => {
      if (!ticking) {
        window.requestAnimationFrame(updateSpace);
        ticking = true;
      }
    };

    // Throttle mouse events
    let mouseThrottle = null;
    const onMouse = (e) => {
      if (mouseThrottle) return;
      mouseThrottle = setTimeout(() => mouseThrottle = null, 50);
      
      mouseX = e.clientX / window.innerWidth;
      mouseY = e.clientY / window.innerHeight;
      if (!ticking) {
        window.requestAnimationFrame(updateSpace);
        ticking = true;
      }
    };

    updateSpace();
    window.addEventListener('scroll', onScroll, { passive: true });
    window.addEventListener('mousemove', onMouse, { passive: true });
    window.addEventListener('resize', updateSpace);
  }

  // MaStar table toggle
  const tableToggle = document.getElementById('mastartable-toggle');
  const tableContainer = document.getElementById('mastartable-container');
  if (tableToggle && tableContainer) {
    tableToggle.addEventListener('click', () => {
      const isOpen = tableToggle.getAttribute('aria-expanded') === 'true';
      const nextState = !isOpen;
      tableToggle.setAttribute('aria-expanded', String(nextState));
      tableContainer.classList.toggle('open', nextState);
    });
  }
});
