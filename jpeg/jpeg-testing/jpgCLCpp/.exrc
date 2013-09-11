if &cp | set nocp | endif
let s:cpo_save=&cpo
set cpo&vim
inoremap <silent> <Plug>(neocomplcache_start_omni_complete) 
inoremap <silent> <Plug>(neocomplcache_start_auto_complete_no_select) 
inoremap <silent> <Plug>(neocomplcache_start_auto_complete) =neocomplcache#mappings#popup_post()
inoremap <silent> <expr> <Plug>(neocomplcache_start_unite_quick_match) unite#sources#neocomplcache#start_quick_match()
inoremap <silent> <expr> <Plug>(neocomplcache_start_unite_complete) unite#sources#neocomplcache#start_complete()
inoremap <silent> <SNR>41_yrrecord =YRRecord3()
inoremap <silent> <expr> <Plug>(neosnippet_start_unite_snippet) unite#sources#snippet#start_complete()
inoremap <silent> <expr> <Plug>(neosnippet_jump) neosnippet#jump_impl()
inoremap <silent> <expr> <Plug>(neosnippet_expand) neosnippet#expand_impl()
inoremap <silent> <expr> <Plug>(neosnippet_jump_or_expand) neosnippet#jump_or_expand_impl()
inoremap <silent> <expr> <Plug>(neosnippet_expand_or_jump) neosnippet#expand_or_jump_impl()
nmap  :GtagsCursor
xmap  :GtagsCursor
omap  :GtagsCursor
snoremap  a<BS>
smap <expr> 	 neosnippet#expandable_or_jumpable() ? "\<Plug>(neosnippet_expand_or_jump)" : "\	"
nmap 	 :Gtags -f %
xmap 	 :Gtags -f %
omap 	 :Gtags -f %
smap  <Plug>(neosnippet_expand_or_jump)
snoremap  a<BS>
nnoremap <silent>  :YRReplace '1', p
nnoremap <silent>  :YRReplace '-1', P
nnoremap ,? ?
nnoremap ,/ /
nnoremap <expr> / ":".v:count1."M/"
nnoremap <expr> ? ":".v:count1."M?"
nmap @ :YRMapsMacro
xmap <silent> K <Plug>(ref-keyword)
nmap <silent> K <Plug>(ref-keyword)
xnoremap <silent> P :YRPaste 'P', 'v'
nnoremap <silent> P :YRPaste 'P'
nmap \rwp <Plug>RestoreWinPosn
xmap \rwp <Plug>RestoreWinPosn
omap \rwp <Plug>RestoreWinPosn
nmap \swp <Plug>SaveWinPosn
xmap \swp <Plug>SaveWinPosn
omap \swp <Plug>SaveWinPosn
nmap \tt <Plug>AM_tt
xmap \tt <Plug>AM_tt
omap \tt <Plug>AM_tt
nmap \tsq <Plug>AM_tsq
xmap \tsq <Plug>AM_tsq
omap \tsq <Plug>AM_tsq
nmap \tsp <Plug>AM_tsp
xmap \tsp <Plug>AM_tsp
omap \tsp <Plug>AM_tsp
nmap \tml <Plug>AM_tml
xmap \tml <Plug>AM_tml
omap \tml <Plug>AM_tml
nmap \tab <Plug>AM_tab
xmap \tab <Plug>AM_tab
omap \tab <Plug>AM_tab
nmap \m= <Plug>AM_m=
xmap \m= <Plug>AM_m=
omap \m= <Plug>AM_m=
nmap \tW@ <Plug>AM_tW@
xmap \tW@ <Plug>AM_tW@
omap \tW@ <Plug>AM_tW@
nmap \t@ <Plug>AM_t@
xmap \t@ <Plug>AM_t@
omap \t@ <Plug>AM_t@
nmap \t~ <Plug>AM_t~
xmap \t~ <Plug>AM_t~
omap \t~ <Plug>AM_t~
nmap \t? <Plug>AM_t?
xmap \t? <Plug>AM_t?
omap \t? <Plug>AM_t?
nmap \w= <Plug>AM_w=
xmap \w= <Plug>AM_w=
omap \w= <Plug>AM_w=
nmap \ts= <Plug>AM_ts=
xmap \ts= <Plug>AM_ts=
omap \ts= <Plug>AM_ts=
nmap \ts< <Plug>AM_ts<
xmap \ts< <Plug>AM_ts<
omap \ts< <Plug>AM_ts<
nmap \ts; <Plug>AM_ts;
xmap \ts; <Plug>AM_ts;
omap \ts; <Plug>AM_ts;
nmap \ts: <Plug>AM_ts:
xmap \ts: <Plug>AM_ts:
omap \ts: <Plug>AM_ts:
nmap \ts, <Plug>AM_ts,
xmap \ts, <Plug>AM_ts,
omap \ts, <Plug>AM_ts,
nmap \t= <Plug>AM_t=
xmap \t= <Plug>AM_t=
omap \t= <Plug>AM_t=
nmap \t< <Plug>AM_t<
xmap \t< <Plug>AM_t<
omap \t< <Plug>AM_t<
nmap \t; <Plug>AM_t;
xmap \t; <Plug>AM_t;
omap \t; <Plug>AM_t;
nmap \t: <Plug>AM_t:
xmap \t: <Plug>AM_t:
omap \t: <Plug>AM_t:
nmap \t, <Plug>AM_t,
xmap \t, <Plug>AM_t,
omap \t, <Plug>AM_t,
nmap \t# <Plug>AM_t#
xmap \t# <Plug>AM_t#
omap \t# <Plug>AM_t#
map \t| <Plug>AM_t|
nmap \T~ <Plug>AM_T~
xmap \T~ <Plug>AM_T~
omap \T~ <Plug>AM_T~
nmap \Tsp <Plug>AM_Tsp
xmap \Tsp <Plug>AM_Tsp
omap \Tsp <Plug>AM_Tsp
nmap \Tab <Plug>AM_Tab
xmap \Tab <Plug>AM_Tab
omap \Tab <Plug>AM_Tab
nmap \TW@ <Plug>AM_TW@
xmap \TW@ <Plug>AM_TW@
omap \TW@ <Plug>AM_TW@
nmap \T@ <Plug>AM_T@
xmap \T@ <Plug>AM_T@
omap \T@ <Plug>AM_T@
nmap \T? <Plug>AM_T?
xmap \T? <Plug>AM_T?
omap \T? <Plug>AM_T?
nmap \T= <Plug>AM_T=
xmap \T= <Plug>AM_T=
omap \T= <Plug>AM_T=
nmap \T< <Plug>AM_T<
xmap \T< <Plug>AM_T<
omap \T< <Plug>AM_T<
nmap \T; <Plug>AM_T;
xmap \T; <Plug>AM_T;
omap \T; <Plug>AM_T;
nmap \T: <Plug>AM_T:
xmap \T: <Plug>AM_T:
omap \T: <Plug>AM_T:
nmap \Ts, <Plug>AM_Ts,
xmap \Ts, <Plug>AM_Ts,
omap \Ts, <Plug>AM_Ts,
nmap \T, <Plug>AM_T,o
xmap \T, <Plug>AM_T,o
omap \T, <Plug>AM_T,o
nmap \T# <Plug>AM_T#
xmap \T# <Plug>AM_T#
omap \T# <Plug>AM_T#
map \T| <Plug>AM_T|
nmap \Htd <Plug>AM_Htd
xmap \Htd <Plug>AM_Htd
omap \Htd <Plug>AM_Htd
nmap \anum <Plug>AM_aunum
xmap \anum <Plug>AM_aunum
omap \anum <Plug>AM_aunum
nmap \aenum <Plug>AM_aenum
xmap \aenum <Plug>AM_aenum
omap \aenum <Plug>AM_aenum
nmap \aunum <Plug>AM_aunum
xmap \aunum <Plug>AM_aunum
omap \aunum <Plug>AM_aunum
nmap \afnc <Plug>AM_afnc
xmap \afnc <Plug>AM_afnc
omap \afnc <Plug>AM_afnc
nmap \adef <Plug>AM_adef
xmap \adef <Plug>AM_adef
omap \adef <Plug>AM_adef
nmap \adec <Plug>AM_adec
xmap \adec <Plug>AM_adec
omap \adec <Plug>AM_adec
nmap \ascom <Plug>AM_ascom
xmap \ascom <Plug>AM_ascom
omap \ascom <Plug>AM_ascom
nmap \aocom <Plug>AM_aocom
xmap \aocom <Plug>AM_aocom
omap \aocom <Plug>AM_aocom
nmap \adcom <Plug>AM_adcom
xmap \adcom <Plug>AM_adcom
omap \adcom <Plug>AM_adcom
nmap \acom <Plug>AM_acom
xmap \acom <Plug>AM_acom
omap \acom <Plug>AM_acom
nmap \abox <Plug>AM_abox
xmap \abox <Plug>AM_abox
omap \abox <Plug>AM_abox
nmap \a( <Plug>AM_a(
xmap \a( <Plug>AM_a(
omap \a( <Plug>AM_a(
nmap \a= <Plug>AM_a=
xmap \a= <Plug>AM_a=
omap \a= <Plug>AM_a=
nmap \a< <Plug>AM_a<
xmap \a< <Plug>AM_a<
omap \a< <Plug>AM_a<
nmap \a, <Plug>AM_a,
xmap \a, <Plug>AM_a,
omap \a, <Plug>AM_a,
nmap \a? <Plug>AM_a?
xmap \a? <Plug>AM_a?
omap \a? <Plug>AM_a?
nmap \r <Plug>(quickrun)
xmap \r <Plug>(quickrun)
omap \r <Plug>(quickrun)
xnoremap <silent> d :YRDeleteRange 'v'
nmap gx <Plug>NetrwBrowseX
nnoremap <silent> gp :YRPaste 'gp'
nnoremap <silent> gP :YRPaste 'gP'
xnoremap <silent> p :YRPaste 'p', 'v'
nnoremap <silent> p :YRPaste 'p'
xnoremap <silent> x :YRDeleteRange 'v'
xnoremap <silent> y :YRYankRange 'v'
snoremap <Left> bi
snoremap <Right> a
snoremap <Del> a<BS>
snoremap <BS> a<BS>
nnoremap <silent> <Plug>NetrwBrowseX :call netrw#NetrwBrowseX(expand("<cWORD>"),0)
nmap <silent> <Plug>RestoreWinPosn :call RestoreWinPosn()
nmap <silent> <Plug>SaveWinPosn :call SaveWinPosn()
nmap <SNR>46_WE <Plug>AlignMapsWrapperEnd
nmap <SNR>46_WS <Plug>AlignMapsWrapperStart
xmap <SNR>46_WS <Plug>AlignMapsWrapperStart
omap <SNR>46_WS <Plug>AlignMapsWrapperStart
nnoremap <silent> <SNR>41_yrrecord :call YRRecord3()
xnoremap <silent> <Plug>(quickrun) :QuickRun -mode v
nnoremap <silent> <Plug>(quickrun) :QuickRun -mode n
nnoremap <silent> <Plug>(quickrun-op) :set operatorfunc=quickrun#operatorg@
xnoremap <silent> <Plug>(ref-keyword) :call ref#K('visual')
nnoremap <silent> <Plug>(ref-keyword) :call ref#K('normal')
xnoremap <silent> <Plug>(neosnippet_register_oneshot_snippet) :call neosnippet#register_oneshot_snippet()
xnoremap <silent> <expr> <Plug>(neosnippet_start_unite_snippet_target) unite#sources#snippet_target#start()
xnoremap <silent> <Plug>(neosnippet_expand_target) :call neosnippet#expand_target()
xnoremap <silent> <Plug>(neosnippet_get_selected_text) :call neosnippet#get_selected_text(visualmode(), 1)
snoremap <silent> <expr> <Plug>(neosnippet_jump) neosnippet#jump_impl()
snoremap <silent> <expr> <Plug>(neosnippet_expand) neosnippet#expand_impl()
snoremap <silent> <expr> <Plug>(neosnippet_jump_or_expand) neosnippet#jump_or_expand_impl()
snoremap <silent> <expr> <Plug>(neosnippet_expand_or_jump) neosnippet#expand_or_jump_impl()
imap <expr> 	 neosnippet#expandable_or_jumpable() ? "\<Plug>(neosnippet_expand_or_jump)" : pumvisible() ? "\" : "\	"
imap  <Plug>(neosnippet_expand_or_jump)
let &cpo=s:cpo_save
unlet s:cpo_save
set autoindent
set backspace=2
set balloondelay=100
set completefunc=neocomplcache#complete#manual_complete
set completeopt=preview,menuone
set expandtab
set fileencodings=ucs-bom,utf-8,default,latin1
set helplang=ja
set hidden
set hlsearch
set ignorecase
set incsearch
set laststatus=2
set ruler
set runtimepath=~/.vim/bundle//.neobundle,~/.vim,~/.vim/bundle/neocomplcache,~/.vim/bundle/neosnippet,~/.vim/bundle/unite.vim,~/.vim/bundle/unite-colorscheme,~/.vim/bundle/unite-font,~/.vim/bundle/vim-ref,~/.vim/bundle/vim-quickrun,~/.vim/bundle/c.vim,~/.vim/bundle/gtags.vim,~/.vim/bundle/Lucius,~/.vim/bundle/grep.vim,~/.vim/bundle/YankRing.vim,~/.vim/bundle/sudo.vim,~/.vim/bundle/smartchr,~/.vim/bundle/ack.vim,~/.vim/bundle/eregex.vim,~/.vim/bundle/Align,~/.vim/bundle/vimdoc-ja,~/software/share/vim/vimfiles,~/software/share/vim/vim74,~/software/share/vim/vimfiles/after,~/.vim/after,~/.vim/bundle/neobundle.vim/
set shiftwidth=2
set showmatch
set smartindent
set smarttab
set softtabstop=2
set textwidth=80
set title
set visualbell
set wildmenu
" vim: set ft=vim :
