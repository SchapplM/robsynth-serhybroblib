% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% fourbarprisDE1
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in fourbarprisDE1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% Ja_rot [3x1]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = fourbarprisDE1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),uint8(0),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_jacobia_rot_sym_varpar: qJ has to be [1x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbarprisDE1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_jacobia_rot_sym_varpar: pkin has to be [3x1] (double)');
Ja_rot=NaN(3,1);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:10:19
	% EndTime: 2020-05-07 09:10:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0; 0; 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:10:19
	% EndTime: 2020-05-07 09:10:20
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (134->19), mult. (75->25), div. (14->6), fcn. (6->2), ass. (0->16)
	t30 = qJ(1) + pkin(3);
	t36 = -pkin(2) + t30;
	t26 = pkin(1) + t36;
	t37 = -pkin(2) - t30;
	t27 = pkin(1) + t37;
	t38 = t26 * t27;
	t25 = pkin(1) - t36;
	t35 = t25 * t38;
	t24 = pkin(1) - t37;
	t33 = t24 * t35;
	t32 = sqrt(-t33);
	t29 = 0.1e1 / t30 ^ 2;
	t28 = 0.1e1 / t30;
	t23 = -pkin(1) ^ 2 + pkin(2) ^ 2 - qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3);
	t22 = 1 / t23 ^ 2;
	t1 = [0; 0; 0.2e1 * ((t29 * t32 / 0.2e1 - t28 / t32 * (-t35 + (t38 + (t26 - t27) * t25) * t24) / 0.4e1) / t23 + (-t30 * t28 - t23 * t29 / 0.2e1) * t32 * t22) / (-t22 * t33 + 0.1e1) * t30;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:10:20
	% EndTime: 2020-05-07 09:10:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN; NaN; NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:10:19
	% EndTime: 2020-05-07 09:10:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (109->17), mult. (54->17), div. (6->4), fcn. (4->2), ass. (0->14)
	t28 = -pkin(3) - qJ(1);
	t26 = -pkin(2) - t28;
	t21 = pkin(1) + t26;
	t27 = -pkin(2) + t28;
	t22 = pkin(1) + t27;
	t29 = t21 * t22;
	t20 = pkin(1) - t26;
	t25 = t20 * t29;
	t19 = pkin(1) - t27;
	t24 = t19 * t25;
	t23 = sqrt(-t24);
	t18 = pkin(1) ^ 2 + pkin(2) ^ 2 - qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3);
	t17 = 1 / t18 ^ 2;
	t1 = [0; 0; (-0.1e1 / t23 * (-t25 + (t29 + (t21 - t22) * t20) * t19) / t18 / 0.2e1 + 0.2e1 * t28 * t23 * t17) / (-t17 * t24 + 0.1e1);];
	Ja_rot = t1;
end