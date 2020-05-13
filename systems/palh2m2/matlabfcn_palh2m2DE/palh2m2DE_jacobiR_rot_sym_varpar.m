% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% JR_rot [9x4]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh2m2DE_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_jacobiR_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m2DE_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_jacobiR_rot_sym_varpar: pkin has to be [5x1] (double)');
JR_rot=NaN(9,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:52
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:51
	% EndTime: 2020-05-03 01:06:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0; t9, 0, 0, 0; 0, 0, 0, 0; -t9, 0, 0, 0; -t8, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:52
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t8 = sin(qJ(2));
	t9 = sin(qJ(1));
	t15 = t9 * t8;
	t11 = cos(qJ(1));
	t14 = t11 * t8;
	t10 = cos(qJ(2));
	t13 = t9 * t10;
	t12 = t11 * t10;
	t1 = [-t13, -t14, 0, 0; t12, -t15, 0, 0; 0, t10, 0, 0; t15, -t12, 0, 0; -t14, -t13, 0, 0; 0, -t8, 0, 0; t11, 0, 0, 0; t9, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:52
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0; t9, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; t9, 0, 0, 0; t8, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:52
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t12 = sin(qJ(3));
	t13 = sin(qJ(1));
	t19 = t13 * t12;
	t14 = cos(qJ(3));
	t18 = t13 * t14;
	t15 = cos(qJ(1));
	t17 = t15 * t12;
	t16 = t15 * t14;
	t1 = [-t18, 0, -t17, 0; t16, 0, -t19, 0; 0, 0, t14, 0; t19, 0, -t16, 0; -t17, 0, -t18, 0; 0, 0, -t12, 0; t15, 0, 0, 0; t13, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:52
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t12 = cos(qJ(1));
	t11 = sin(qJ(1));
	t1 = [-t11, 0, 0, 0; t12, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; t12, 0, 0, 0; t11, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:52
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->4), mult. (16->4), div. (0->0), fcn. (32->4), ass. (0->7)
	t35 = cos(qJ(1));
	t34 = cos(qJ(4));
	t33 = sin(qJ(1));
	t32 = sin(qJ(4));
	t30 = -t35 * t32 - t33 * t34;
	t29 = t33 * t32 - t35 * t34;
	t1 = [t30, 0, 0, t30; -t29, 0, 0, -t29; 0, 0, 0, 0; t29, 0, 0, t29; t30, 0, 0, t30; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
end