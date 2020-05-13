% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh2m1OL_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m1OL_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_jacobiR_rot_sym_varpar: pkin has to be [6x1] (double)');
JR_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:40
	% EndTime: 2020-05-03 00:53:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->9), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t8 = sin(qJ(2));
	t9 = sin(qJ(1));
	t15 = t9 * t8;
	t11 = cos(qJ(1));
	t14 = t11 * t8;
	t10 = cos(qJ(2));
	t13 = t9 * t10;
	t12 = t11 * t10;
	t1 = [-t13, -t14, 0, 0, 0; t12, -t15, 0, 0, 0; 0, -t10, 0, 0, 0; t15, -t12, 0, 0, 0; -t14, -t13, 0, 0, 0; 0, t8, 0, 0, 0; -t11, 0, 0, 0, 0; -t9, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->15), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t21 = qJ(2) + qJ(3);
	t19 = sin(t21);
	t22 = sin(qJ(1));
	t27 = t22 * t19;
	t20 = cos(t21);
	t26 = t22 * t20;
	t23 = cos(qJ(1));
	t25 = t23 * t19;
	t24 = t23 * t20;
	t1 = [-t26, -t25, -t25, 0, 0; t24, -t27, -t27, 0, 0; 0, -t20, -t20, 0, 0; t27, -t24, -t24, 0, 0; -t25, -t26, -t26, 0, 0; 0, t19, t19, 0, 0; -t23, 0, 0, 0, 0; -t22, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (63->20), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t28 = qJ(2) + qJ(3) + qJ(4);
	t26 = sin(t28);
	t29 = sin(qJ(1));
	t34 = t29 * t26;
	t27 = cos(t28);
	t33 = t29 * t27;
	t30 = cos(qJ(1));
	t32 = t30 * t26;
	t31 = t30 * t27;
	t1 = [-t33, -t32, -t32, -t32, 0; t31, -t34, -t34, -t34, 0; 0, -t27, -t27, -t27, 0; t34, -t31, -t31, -t31, 0; -t32, -t33, -t33, -t33, 0; 0, t26, t26, t26, 0; -t30, 0, 0, 0, 0; -t29, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (99->20), mult. (64->20), div. (0->0), fcn. (111->6), ass. (0->24)
	t104 = qJ(2) + qJ(3) + qJ(4);
	t103 = cos(t104);
	t107 = cos(qJ(5));
	t115 = t103 * t107;
	t105 = sin(qJ(5));
	t106 = sin(qJ(1));
	t114 = t106 * t105;
	t113 = t106 * t107;
	t108 = cos(qJ(1));
	t112 = t108 * t105;
	t111 = t108 * t107;
	t102 = sin(t104);
	t110 = t102 * t113;
	t109 = t102 * t111;
	t101 = t108 * t103;
	t100 = t103 * t105;
	t99 = t106 * t103;
	t98 = t102 * t112;
	t97 = t102 * t114;
	t96 = t103 * t111 - t114;
	t95 = -t103 * t112 - t113;
	t94 = -t103 * t113 - t112;
	t93 = t103 * t114 - t111;
	t1 = [t94, -t109, -t109, -t109, t95; t96, -t110, -t110, -t110, -t93; 0, -t115, -t115, -t115, t102 * t105; t93, t98, t98, t98, -t96; t95, t97, t97, t97, t94; 0, t100, t100, t100, t102 * t107; -t106 * t102, t101, t101, t101, 0; t108 * t102, t99, t99, t99, 0; 0, -t102, -t102, -t102, 0;];
	JR_rot = t1;
end