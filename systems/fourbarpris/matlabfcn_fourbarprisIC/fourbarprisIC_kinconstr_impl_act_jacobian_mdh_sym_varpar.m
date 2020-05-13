% Jacobian of implicit kinematic constraints of
% fourbarprisIC
% with respect to active joint coordinates
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
% Phi_a [(no of constraints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:59
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi_a = fourbarprisIC_kinconstr_impl_act_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisIC_kinconstr_impl_act_jacobian_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisIC_kinconstr_impl_act_jacobian_mdh_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_active_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:59:30
% EndTime: 2020-05-07 09:59:30
% DurationCPUTime: 0.01s
% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
t1 = [-cos(qJ(1)); -sin(qJ(1)); 0;];
Phi_a = t1;
