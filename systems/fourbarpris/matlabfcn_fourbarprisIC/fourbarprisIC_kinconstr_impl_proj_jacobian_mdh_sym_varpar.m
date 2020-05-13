% Jacobian of implicit kinematic constraints of
% fourbarprisIC
% projection from active to passive joints coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% B21 [(no of passive joints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:59
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function B21 = fourbarprisIC_kinconstr_impl_proj_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisIC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisIC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_projection_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:59:30
% EndTime: 2020-05-07 09:59:30
% DurationCPUTime: 0.02s
% Computational Cost: add. (10->7), mult. (11->11), div. (7->4), fcn. (18->4), ass. (0->8)
t21 = sin(qJ(3));
t22 = sin(qJ(1));
t23 = cos(qJ(3));
t24 = cos(qJ(1));
t27 = 0.1e1 / pkin(2) / (t21 * t24 - t22 * t23);
t26 = pkin(3) + qJ(2);
t20 = 0.1e1 / t26;
t1 = [(-t23 * t21 - t24 * t22) * t20 / (-t23 ^ 2 + t24 ^ 2); -t27; ((-t21 * t22 - t23 * t24) * pkin(2) + t26) * t20 * t27;];
B21 = t1;
