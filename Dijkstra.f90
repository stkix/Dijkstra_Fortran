program Dijkstra_Main
    implicit none
    integer :: num_vertices, num_edges, i, start_vertex, end_vertex, count
    !num_verticesは点の数、num_edgesは辺の数、start_vertexは始点、end_vertexは終点
    real(8), allocatable :: x(:), y(:), weight(:,:)
    !weightは(重み)（グラフ理論）
    real(8), allocatable :: dist(:)
    !distは最短距離を収納するために使われる
    integer, allocatable :: previous(:)
    !各頂点の最短経路における直前の頂点のインデックス（索引）を収納する
  
    ! 頂点数を数える
    open(unit=10, file='inp.txt', status="old", action="read")
    num_vertices = 0
    do
      read(10, *, iostat=count)
      if (count == 0) then
        num_vertices = num_vertices + 1
      else
        exit
      end if
    end do
  
   !ファイルの読み取りの初期化
    rewind(unit=10)
    
    ! 配列を割り当てる
    allocate(x(1:num_vertices))
    allocate(y(1:num_vertices))
    allocate(weight(1:num_vertices,1:num_vertices))
    allocate(dist(1:num_vertices))
    allocate(previous(1:num_vertices))
  
    ! xy座標を読み取る
    do i = 1, num_vertices
      read(10, *) x(i), y(i)
    end do
  
    close(unit=10)
  
    ! 辺の数を入力する
    write(*,*) "頂点数: ", num_vertices
    write(*,'(100("-"))')
    write(*,*) "辺数を入力してください:"
    read(*,*) num_edges
    write(*,'(100("-"))')
  
    ! 辺の数が正しいかチェックする
    if (num_edges < 1) then
      write(*,*) "エラー: 辺数は1以上にしなければならない"
      write(*,'(100("-"))')
      stop
    end if
  
    ! 辺の情報を入力する
    call read_edges(num_edges, num_vertices, x, y, weight)
  
    ! 始点と終点を入力する
    do
      write(*,*) "始点を入力してください (1 から ", num_vertices, "の範囲内):"
      read(*,*) start_vertex
      write(*,'(100("-"))')
      if (start_vertex >= 1 .and. start_vertex <= num_vertices) exit
      write(*,*) "エラー: 無効な始点番号が指定されました"
      write(*,'(100("-"))')
    end do
    do
      write(*,*) "終点を入力してください (1 から ", num_vertices, "の範囲内):"
      read(*,*) end_vertex
      write(*,'(100("-"))')
      if (end_vertex >= 1 .and. end_vertex <= num_vertices) exit
      write(*,*) "エラー: 無効な終点番号が指定されました"
      write(*,'(100("-"))')
    end do
  
    ! ダイクストラ法で最短経路を計算する
    call dijkstra_algorithm(num_vertices, start_vertex, weight, dist, previous)
  
    ! 結果を表示する
    do i = 1, num_vertices
      write(*,*) "頂点 ", i, " までの最短距離: ", dist(i)
      write(*,'(100("-"))')
    end do
  
    ! 最短経路を表示する
    do i = 1, num_vertices
      write(*,*) "頂点 ", i, " への最短経路:"
      write(*,*) write_shortest_path(i, previous)
      write(*,'(100("-"))')
    end do
  
  
    !終点の情報を表示する
      write(*,*) "終点 ", end_vertex, " までの最短距離: ", dist(end_vertex)
      write(*,'(100("-"))')
  
      write(*,*) "終点 ", end_vertex, " への最短経路:"
      write(*,*) write_shortest_path(end_vertex, previous)
      write(*,'(100("-"))')
  
    ! メモリを解放する
    deallocate(x, y, weight, dist, previous)
  
  
  
  
  contains
  
  
  
  
    ! 辺の情報を読み取るサブルーチン
    subroutine read_edges(num_edges, num_vertices, x, y, weight)
      implicit none
      integer, intent(in) :: num_edges, num_vertices
      real(8), intent(in) :: x(:), y(:)
      real(8), allocatable, intent(out) :: weight(:,:)
      integer :: i, start_vertex, end_vertex
      real(8) :: infinity
  
      ! 重み配列を大きな値 (無限大を表す) で初期化する。
      infinity = 1.0e30
      allocate(weight(1:num_vertices, 1:num_vertices))
      weight = infinity
  
      write(*,*) "各辺の始点と終点を入力してください:"
      write(*,'(100("-"))')
      do i = 1, num_edges
        write(*,*) "辺 ", i
        read(*,*) start_vertex, end_vertex
        if (start_vertex < 1 .or. start_vertex > num_vertices .or. end_vertex < 1 .or. end_vertex > num_vertices) then
          write(*,*) "エラー: 無効な頂点番号が指定されました"
          write(*,'(100("-"))')
          stop
        end if
        weight(start_vertex, end_vertex) = sqrt((x(start_vertex) - x(end_vertex))**2 + (y(start_vertex) - y(end_vertex))**2)
        weight(end_vertex, start_vertex) = weight(start_vertex, end_vertex) ! 無向グラフの場合、両方の方向に重みを設定する
        !有向グラフにおいては辺1や辺2などは一つの方向性しか持たず、一方通行しかできないが、無向グラフにおいては頂点と頂点との相互方向性（二つでお互い逆の方向）があり、相互通行を行うことができる。
        !この問題では点と点の距離を求めたいので求めるのはスカラー量である。そのためこの場合は無向グラフが問題に適していると考えられる。
        write(*,*) "辺 ", i, " の重み（距離）は",weight(end_vertex, start_vertex),"である"
        write(*,'(100("-"))')
      end do
    end subroutine read_edges
  
    ! ダイクストラ法のサブルーチン
    subroutine dijkstra_algorithm(num_vertices, start_vertex, weight, dist, previous)
      implicit none
      integer, intent(in) :: num_vertices, start_vertex
      real(8), intent(in) :: weight(:,:)
      real(8), allocatable, intent(out) :: dist(:)
      integer, allocatable, intent(out) :: previous(:)
      logical, allocatable :: visited(:)!点を訪問しているかいないかの時に使われる
      integer :: i, j, current_vertex
      !current_vertexは現在の点を示す
      real(8) :: infinity!距離の初期値を設定するのに使われる
  
      infinity = 1.0e30
  
      allocate(visited(1:num_vertices))
      allocate(dist(1:num_vertices))
      allocate(previous(1:num_vertices))
  
      visited = .false.
      dist = infinity
      previous = -1  ! 'previous' 配列を無効な状態に初期化する
  
      dist(start_vertex) = 0.0  !始点の最短距離を0とする。
  
      do i = 1, num_vertices
        current_vertex = 0
  
        !すべて未確認の点のうち距離が一番短い点を確定する。（現在加算された重みが一番小さい点）
        do j = 1, num_vertices
          if (.not. visited(j) .and. (current_vertex == 0 .or. dist(j) < dist(current_vertex))) then
            current_vertex = j
          end if
        end do
  
        visited(current_vertex) = .true.
  
        !先ほど確定した点から未確定の点までの距離を計算し、もしそれが最短距離の場合は更新する。
        do j = 1, num_vertices
          if (.not. visited(j) .and. dist(current_vertex) + weight(current_vertex, j) < dist(j)) then
            dist(j) = dist(current_vertex) + weight(current_vertex, j)
            previous(j) = current_vertex
          end if
        end do
      end do
  
      deallocate(visited)
  
    end subroutine dijkstra_algorithm
  
  ! 最短経路を表示する関数副プログラム
  ! 目標の点から始点に向かって逆順に最短経路をたどり、それを配列として返す。
    function write_shortest_path(vertex, previous) result(path)
      implicit none
      integer, intent(in) :: vertex ! 最短経路を求める終点のインデックス（索引）
      integer, intent(in) :: previous(:) ! 各頂点の最短経路における直前の頂点のインデックス（索引）を収納する
      integer :: next_vertex, i, num_vertices
      integer, allocatable :: path(:)
    
      num_vertices = size(previous)
      allocate(path(num_vertices))
      next_vertex = vertex
      path(1) = next_vertex
      i = 2
    
      do while (previous(next_vertex) /= -1)
        next_vertex = previous(next_vertex)
        path(i) = next_vertex
        i = i + 1
      end do
    
      path = path(1:i-1)
    
      ! pathを逆転する
      call reverse_array(path)
    
    end function write_shortest_path
    
    ! 配列を逆転するサブルーチン
    subroutine reverse_array(arr)
      implicit none
      integer, intent(inout) :: arr(:)
      integer :: i, temp
    
      do i = 1, size(arr) / 2
        temp = arr(i)
        arr(i) = arr(size(arr) - i + 1)
        arr(size(arr) - i + 1) = temp
      end do
    
    end subroutine reverse_array
  
  end program Dijkstra_Main